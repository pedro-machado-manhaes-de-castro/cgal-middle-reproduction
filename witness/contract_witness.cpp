// Table: exhaustive six-point ordering checks in README.md.
#include <CGAL/Simple_cartesian.h>
#include <CGAL/hilbert_sort.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

using K = CGAL::Simple_cartesian<double>;
using Point = K::Point_2;

std::vector<Point> make_points(const std::array<int, 6>& numerators)
{
    std::vector<Point> points;
    for(int value : numerators)
        points.emplace_back(value / 64.0, value / 64.0);
    return points;
}

bool ordered(const std::vector<Point>& points)
{
    return std::is_sorted(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a.x() < b.x();
    });
}

int main()
{
    const std::array<int, 6> example = {0, 5, 9, 7, 11, 64};
    auto sample = make_points(example);
    CGAL::Hilbert_sort_2<K, CGAL::Hilbert_sort_middle_policy>()(sample.begin(), sample.end());
    std::cout << "{\"example_output_numerators\":[";
    for(std::size_t i = 0; i < sample.size(); ++i)
        std::cout << (i ? "," : "") << static_cast<int>(sample[i].x() * 64);
    std::cout << "]";

    std::array<int, 6> values = {0, 5, 7, 9, 11, 64};
    unsigned permutations = 0, direct_failures = 0, function_failures = 0;
    do {
        auto direct = make_points(values);
        auto function = direct;
        CGAL::Hilbert_sort_2<K, CGAL::Hilbert_sort_middle_policy>()(direct.begin(), direct.end());
        CGAL::hilbert_sort(function.begin(), function.end(), CGAL::Hilbert_sort_middle_policy());
        direct_failures += !ordered(direct);
        function_failures += !ordered(function);
        ++permutations;
    } while(std::next_permutation(values.begin(), values.end()));
    std::cout << ",\"permutations\":" << permutations
              << ",\"direct_failures\":" << direct_failures
              << ",\"function_failures\":" << function_failures << "}\n";
}
