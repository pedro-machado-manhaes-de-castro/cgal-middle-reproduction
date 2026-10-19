#include <CGAL/Simple_cartesian.h>
#include <CGAL/hilbert_sort.h>

#include <iostream>
#include <vector>

int main()
{
    using Kernel = CGAL::Simple_cartesian<double>;
    using Point = Kernel::Point_2;

    std::vector<Point> points;
    for (int numerator : {0, 5, 9, 7, 11, 64})
        points.emplace_back(numerator / 64.0, numerator / 64.0);

    CGAL::Hilbert_sort_2<Kernel, CGAL::Hilbert_sort_middle_policy> sort;
    sort(points.begin(), points.end());

    for (const Point& point : points)
        std::cout << point.x() * 64 << ' ';
    std::cout << '\n';
}
