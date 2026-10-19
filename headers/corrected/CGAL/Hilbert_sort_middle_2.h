// Copyright (c) 2011  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     :  Olivier Devillers

#ifndef CGAL_HILBERT_SORT_MIDDLE_2_H
#define CGAL_HILBERT_SORT_MIDDLE_2_H

#include <CGAL/config.h>
#include <functional>
#include <cstddef>
#include <CGAL/Hilbert_sort_middle_base.h>
#include <CGAL/number_utils.h>

namespace CGAL {

namespace internal {
template <class K, int x, bool up> struct Fixed_hilbert_cmp_2;

template <class K, int x>
struct Fixed_hilbert_cmp_2<K,x,true>
    : public CGAL::cpp98::binary_function<typename K::Point_2,
    typename K::Point_2, bool>
{
  typedef typename K::Point_2 Point;
  K k;
  double value;
  Fixed_hilbert_cmp_2 (double v, const K &_k) : k(_k),value(v) {}
  bool operator() (const Point &p) const
  {
    return ! Fixed_hilbert_cmp_2<K,x,false> (value, k) (p);
  }
};

template <class K>
struct Fixed_hilbert_cmp_2<K,0,false>
    : public CGAL::cpp98::binary_function<typename K::Point_2,
    typename K::Point_2, bool>
{
  typedef typename K::Point_2 Point;
  K k;
  double value;
  Fixed_hilbert_cmp_2 (double v, const K &_k) : k(_k),value(v) {}
  bool operator() (const Point &p) const
  {
    return to_double(k.compute_x_2_object()(p)) < value;
  }
};

template <class K>
struct Fixed_hilbert_cmp_2<K,1,false>
    : public CGAL::cpp98::binary_function<typename K::Point_2,
    typename K::Point_2, bool>
{
  typedef typename K::Point_2 Point;
  K k;
  double value;
  Fixed_hilbert_cmp_2 (double v, const K &_k) : k(_k),value(v) {}
  bool operator() (const Point &p) const
  {
    return to_double(k.compute_y_2_object()(p)) < value;
  }
};

} // namespace internal

template <class K>
class Hilbert_sort_middle_2
{
public:
  typedef K Kernel;
  typedef typename Kernel::Point_2 Point;

private:
  Kernel _k;
  std::ptrdiff_t _limit;

  template <int x, bool up>
  struct Cmp
    : public internal::Fixed_hilbert_cmp_2<Kernel,x,up>
  {
    Cmp (double v, const Kernel &k) : internal::Fixed_hilbert_cmp_2<Kernel,x,up> (v, k) {}
  };

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  void skip_empty_cells(RandomAccessIterator begin, RandomAccessIterator end,
                        double xmin, double ymin, double xmax, double ymax,
                        double pxmin, double pymin, double pxmax, double pymax) const
  {
    const int y = (x + 1) % 2;
    const double xmed = (xmin + xmax) / 2;
    const double ymed = (ymin + ymax) / 2;
    if ((pxmin < xmed && !(pxmax < xmed)) ||
        (pymin < ymed && !(pymax < ymed))) {
      sort<x, upx, upy>(begin, end, xmin, ymin, xmax, ymax);
      return;
    }
    const bool progress_x = (xmin < xmed && xmed < xmax) || (xmax < xmed && xmed < xmin);
    const bool progress_y = (ymin < ymed && ymed < ymax) || (ymax < ymed && ymed < ymin);
    if (!progress_x && !progress_y) return;

    const bool first_x = (pxmin < xmed) != upx;
    const bool first_y = (pymin < ymed) != upy;
    if (first_x) {
      if (first_y)
        skip_empty_cells<y, upy, upx>(begin, end, ymin, xmin, ymed, xmed,
                                     pymin, pxmin, pymax, pxmax);
      else
        skip_empty_cells<x, upx, upy>(begin, end, xmin, ymed, xmed, ymax,
                                     pxmin, pymin, pxmax, pymax);
    } else {
      if (!first_y)
        skip_empty_cells<x, upx, upy>(begin, end, xmed, ymed, xmax, ymax,
                                     pxmin, pymin, pxmax, pymax);
      else
        skip_empty_cells<y, !upy, !upx>(begin, end, ymed, xmax, ymin, xmed,
                                       pymin, pxmin, pymax, pxmax);
    }
  }

  template <class RandomAccessIterator>
  bool compute_bounds(RandomAccessIterator begin, RandomAccessIterator end,
                      double& xmin, double& ymin, double& xmax, double& ymax) const
  {
    xmin = xmax = to_double(_k.compute_x_2_object()(*begin));
    ymin = ymax = to_double(_k.compute_y_2_object()(*begin));
    RandomAccessIterator it = begin + 1;
    for(; it < end; ++it)
      if (to_double(_k.compute_x_2_object()(*it)) != xmin ||
          to_double(_k.compute_y_2_object()(*it)) != ymin)
        break;
    if (it == end) return false;
    for(; it < end; ++it) {
      const double px = to_double(_k.compute_x_2_object()(*it));
      const double py = to_double(_k.compute_y_2_object()(*it));
      if (px < xmin) xmin = px;
      if (px > xmax) xmax = px;
      if (py < ymin) ymin = py;
      if (py > ymax) ymax = py;
    }
    return true;
  }

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  void refine_same_range(RandomAccessIterator begin, RandomAccessIterator end,
                         double xmin, double ymin, double xmax, double ymax) const
  {
    double pxmin, pymin, pxmax, pymax;
    if (!compute_bounds(begin, end, pxmin, pymin, pxmax, pymax)) return;
    if (x == 0)
      skip_empty_cells<x, upx, upy>(begin, end, xmin, ymin, xmax, ymax,
                                   pxmin, pymin, pxmax, pymax);
    else
      skip_empty_cells<x, upx, upy>(begin, end, xmin, ymin, xmax, ymax,
                                   pymin, pxmin, pymax, pxmax);
  }

public:
  Hilbert_sort_middle_2 (const Kernel &k, std::ptrdiff_t limit = 1)
    : _k(k), _limit (limit)
  {}

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  void sort (RandomAccessIterator begin, RandomAccessIterator end,
             double xmin, double ymin, double xmax, double ymax) const
  {
    const int y = (x + 1) % 2;
    if (end - begin <= _limit) return;

    double xmed= (xmin+xmax)/2;
    double ymed= (ymin+ymax)/2;

    RandomAccessIterator m0 = begin, m4 = end;

    RandomAccessIterator m2 = internal::fixed_hilbert_split (m0, m4, Cmp< x,  upx> (xmed,_k));
    RandomAccessIterator m1 = internal::fixed_hilbert_split (m0, m2, Cmp< y,  upy> (ymed,_k));
    RandomAccessIterator m3 = internal::fixed_hilbert_split (m2, m4, Cmp< y, !upy> (ymed,_k));

    if (m1 != m4)
      sort<y, upy, upx>(m0, m1, ymin, xmin, ymed, xmed);
    else
      refine_same_range<y, upy, upx>(m0, m1, ymin, xmin, ymed, xmed);
    if (m1 != m0 || m2 != m4)
      sort<x, upx, upy>(m1, m2, xmin, ymed, xmed, ymax);
    else
      refine_same_range<x, upx, upy>(m1, m2, xmin, ymed, xmed, ymax);
    if (m2 != m0 || m3 != m4)
      sort<x, upx, upy>(m2, m3, xmed, ymed, xmax, ymax);
    else
      refine_same_range<x, upx, upy>(m2, m3, xmed, ymed, xmax, ymax);
    if (m3 != m0)
      sort<y, !upy, !upx>(m3, m4, ymed, xmax, ymin, xmed);
    else
      refine_same_range<y, !upy, !upx>(m3, m4, ymed, xmax, ymin, xmed);
  }

  template <class RandomAccessIterator>
  void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
  {
    //Bbox_2 box=bbox_2(begin, end); BUG: WE NEED TO FIX THIS
    double xmin, ymin, xmax, ymax;
    if (!compute_bounds(begin, end, xmin, ymin, xmax, ymax)) return;
    sort <0, false, false> (begin, end, xmin, ymin, xmax, ymax);
  }
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MIDDLE_2_H
