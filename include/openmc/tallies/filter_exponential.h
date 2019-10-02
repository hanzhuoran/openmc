#ifndef OPENMC_TALLIES_FILTER_EXPONENTIAL_H
#define OPENMC_TALLIES_FILTER_EXPONENTIAL_H

#include <string>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Gives exponential moments of a particle's position
//==============================================================================

class ExponentialFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~ExponentialFilter() = default;

  //----------------------------------------------------------------------------
  // Methods
  std::string type() const override {return "exponential";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int order() const { return order_; }
  virtual void set_order(int order);

  double exponent() const { return exponent_; }
  virtual void set_exponent(double exponent);

  double x() const { return x_; }
  void set_x(double x) { x_ = x; }

  double y() const { return y_; }
  void set_y(double y) { y_ = y; }

  double r() const { return r_; }
  void set_r(double r) { r_ = r; }

  //----------------------------------------------------------------------------
  // Data members

protected:
  //! Cartesian x coordinate for the origin of this expansion.
  double x_;

  //! Cartesian y coordinate for the origin of this expansion.
  double y_;

  //! Maximum radius from the origin covered by this expansion.
  double r_;

  int order_;

  double exponent_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_ZERNIKE_H
