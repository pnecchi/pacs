//----------------------------------------------------------------------
// Description: Vector norms
// Author:      Pierpaolo Necchi
// Email:       pierpaolo.necchi@gmail.com
// Date:        lun 11 apr 2016 14:49:19 CEST
//----------------------------------------------------------------------

#ifndef NORMS_H
#define NORMS_H

#include <vector>

double normSup(std::vector<double> const &v);
double normL2(std::vector<double> const &v, double const h);
double normH1(std::vector<double> const &v, double const h);

#endif /* end of include guard: NORMS_H */
