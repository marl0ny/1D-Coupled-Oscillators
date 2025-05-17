#include <vector>

#ifndef _ORTHOGONAL_TRANSFORMS_
#define _ORTHOGONAL_TRANSFORMS_

void make_dst(std::vector<double> &dst, std::vector<double> &idst, int n);

void make_dsct(std::vector<double> &dsct, std::vector<double> &idsct, int n);

std::vector<double> apply_transform(
    const std::vector<double> &transform, const std::vector<double> &x);

#endif