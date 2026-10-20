// Table: Delaunay insertion paths and phase times in analysis/paired_ratios.csv.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/spatial_sort.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;
using Triangulation = CGAL::Delaunay_triangulation_2<Kernel>;
using Clock = std::chrono::steady_clock;

std::uint64_t mix(std::uint64_t value)
{
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

std::uint64_t point_hash(const Point& point)
{
    double x = point.x(), y = point.y();
    std::uint64_t a, b;
    std::memcpy(&a, &x, sizeof(a));
    std::memcpy(&b, &y, sizeof(b));
    return mix(a) ^ mix(b + 0x9e3779b97f4a7c15ULL);
}

std::uint64_t multiset_hash(const std::vector<Point>& points)
{
    std::uint64_t hash = 0;
    for(const auto& point : points)
        hash += point_hash(point);
    return hash;
}

std::uint64_t order_hash(const std::vector<Point>& points)
{
    std::uint64_t hash = 1469598103934665603ULL;
    for(const auto& point : points)
        hash = (hash ^ point_hash(point)) * 1099511628211ULL;
    return hash;
}

long double path_length(const std::vector<Point>& points)
{
    long double length = 0;
    for(std::size_t i = 1; i < points.size(); ++i) {
        const long double dx = static_cast<long double>(points[i].x()) - points[i - 1].x();
        const long double dy = static_cast<long double>(points[i].y()) - points[i - 1].y();
        length += std::hypot(dx, dy);
    }
    return length;
}

std::size_t unique_count(std::vector<Point> points)
{
    std::sort(points.begin(), points.end(), Kernel().less_xy_2_object());
    return std::distance(points.begin(), std::unique(points.begin(), points.end()));
}

std::vector<Point> generate(const std::string& family, std::size_t n, unsigned seed)
{
    std::mt19937_64 random(seed);
    const auto unit = [&]() { return (random() >> 11) * 0x1.0p-53; };
    std::vector<Point> points;
    points.reserve(n);
    const auto side = static_cast<std::size_t>(std::ceil(std::sqrt(n)));
    for(std::size_t i = 0; i < n; ++i) {
        double x = unit(), y = unit();
        unit();
        if(family == "grid") {
            x = double(i % side) / side;
            y = double(i / side) / side;
        } else if(family == "collinear") {
            y = x;
        } else if(family == "clustered") {
            const double center = (i % 8) / 8.0;
            x = center + x * 0x1p-10;
            y = center + y * 0x1p-10;
        } else if(family == "concentrated") {
            x *= 0x1p-40;
            y *= 0x1p-40;
            if(i == 0)
                x = y = 1;
        } else if(family == "duplicates") {
            x = y = 0.25;
            if(i == 0)
                x = y = 0;
            if(i == 1)
                x = y = 1;
        } else if(family == "exponential") {
            x = std::ldexp(1.0 + x * 0.25, -static_cast<int>(i % 1000));
            y = x;
        } else if(family == "graded") {
            x = std::ldexp(x, -static_cast<int>(i % 40));
            y = std::ldexp(y, -static_cast<int>(i % 40));
        } else if(family != "uniform") {
            throw std::runtime_error("Unknown synthetic family");
        }
        points.emplace_back(x, y);
    }
    std::shuffle(points.begin(), points.end(), random);
    return points;
}

std::vector<Point> load_points(const std::string& file, std::size_t n, unsigned seed)
{
    std::ifstream stream(file, std::ios::binary | std::ios::ate);
    if(!stream || stream.tellg() < 0 || std::uint64_t(stream.tellg()) % 16 != 0)
        throw std::runtime_error("Invalid coordinate file");
    const std::size_t count = std::uint64_t(stream.tellg()) / 16;
    if(n > count)
        throw std::runtime_error("Requested sample exceeds dataset");
    stream.seekg(0);
    std::vector<Point> points;
    points.reserve(count);
    for(std::size_t i = 0; i < count; ++i) {
        double xy[2];
        stream.read(reinterpret_cast<char*>(xy), sizeof(xy));
        if(!stream || !std::isfinite(xy[0]) || !std::isfinite(xy[1]))
            throw std::runtime_error("Invalid coordinate");
        points.emplace_back(xy[0], xy[1]);
    }
    std::mt19937_64 random(seed);
    std::shuffle(points.begin(), points.end(), random);
    points.resize(n);
    return points;
}

void sort_points(std::vector<Point>& points, const std::string& policy)
{
    if(policy == "middle")
        CGAL::spatial_sort(points.begin(), points.end(), CGAL::Hilbert_sort_middle_policy());
    else if(policy == "median")
        CGAL::spatial_sort(points.begin(), points.end(), CGAL::Hilbert_sort_median_policy());
    else
        throw std::runtime_error("Unknown sorting policy");
}

int main(int argc, char** argv)
{
    try {
        if(argc != 7 && argc != 8)
            throw std::runtime_error("Usage: bench policy source n seed repeats verify [order.csv]");
        const std::string policy = argv[1], source = argv[2];
        const std::size_t n = std::stoull(argv[3]);
        const unsigned seed = static_cast<unsigned>(std::stoul(argv[4]));
        const unsigned repeats = static_cast<unsigned>(std::stoul(argv[5]));
        const bool verify = std::stoi(argv[6]) != 0;
        if(n < 2 || repeats < 1)
            throw std::runtime_error("At least two points and one repetition are required");
        const auto input = source.rfind("synthetic:", 0) == 0
            ? generate(source.substr(10), n, seed) : load_points(source, n, seed);
        const auto input_hash = multiset_hash(input);
        const auto input_order = order_hash(input);
        const auto input_length = path_length(input);
        const auto expected_vertices = unique_count(input);
        auto warm = generate("uniform", 4096, 987654);
        sort_points(warm, policy);
        {
            Triangulation triangulation;
            Triangulation::Face_handle hint;
            for(const auto& point : warm)
                hint = triangulation.insert(point, hint)->face();
        }
        for(unsigned trial = 0; trial < repeats; ++trial) {
            auto points = input;
            Triangulation triangulation;
            Triangulation::Face_handle hint;
#ifdef BENCH_TRACE_INSERTION
            std::vector<Point> inserted;
            inserted.reserve(n);
#endif
            const auto cpu_start = std::clock();
            const auto start = Clock::now();
            sort_points(points, policy);
            const auto split = Clock::now();
            const auto cpu_split = std::clock();
            for(const auto& point : points) {
#ifdef BENCH_TRACE_INSERTION
                inserted.push_back(point);
#endif
                hint = triangulation.insert(point, hint)->face();
            }
            const auto finish = Clock::now();
            const auto cpu_finish = std::clock();
            const double sorting = std::chrono::duration<double>(split - start).count();
            const double insertion = std::chrono::duration<double>(finish - split).count();
            const double total = std::chrono::duration<double>(finish - start).count();
            const auto output_hash = multiset_hash(points);
            const auto output_order = order_hash(points);
            const auto output_length = path_length(points);
            const bool valid = !verify || triangulation.is_valid();
            if(output_hash != input_hash || !valid || !std::isfinite(output_length)
               || triangulation.number_of_vertices() != expected_vertices)
                throw std::runtime_error("Construction validation failed");
#ifdef BENCH_TRACE_INSERTION
            if(inserted != points)
                throw std::runtime_error("Insertion trace differs from the measured sequence");
#endif
            if(argc == 8 && trial == 0) {
                std::ofstream dump(argv[7]);
                dump << std::setprecision(17) << "x,y\n";
                for(const auto& point : points)
                    dump << point.x() << ',' << point.y() << '\n';
                if(!dump)
                    throw std::runtime_error("Could not write insertion sequence");
            }
            std::cout << std::setprecision(20)
                      << "{\"variant\":\"" << BENCH_VARIANT << "\",\"policy\":\"" << policy
                      << "\",\"n\":" << n << ",\"seed\":" << seed << ",\"trial\":" << trial
                      << ",\"sort_seconds\":" << sorting << ",\"insert_seconds\":" << insertion
                      << ",\"total_seconds\":" << total
                      << ",\"sort_cpu_seconds\":" << double(cpu_split - cpu_start) / CLOCKS_PER_SEC
                      << ",\"insert_cpu_seconds\":" << double(cpu_finish - cpu_split) / CLOCKS_PER_SEC
                      << ",\"total_cpu_seconds\":" << double(cpu_finish - cpu_start) / CLOCKS_PER_SEC
                      << ",\"input_path\":" << input_length << ",\"insertion_path\":" << output_length
                      << ",\"input_hash\":\"" << input_hash << "\",\"output_hash\":\"" << output_hash
                      << "\",\"input_order_hash\":\"" << input_order << "\",\"insertion_order_hash\":\"" << output_order
                      << "\",\"vertices\":" << triangulation.number_of_vertices()
                      << ",\"faces\":" << triangulation.number_of_faces()
                      << ",\"expected_vertices\":" << expected_vertices
                      << ",\"dimension\":" << triangulation.dimension()
                      << ",\"validity_checked\":" << (verify ? "true" : "false")
                      << ",\"valid\":" << (valid ? "true" : "false") << "}\n" << std::flush;
        }
    } catch(const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
