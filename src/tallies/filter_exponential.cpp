#include "openmc/tallies/filter_exponential.h"

#include <cmath>
#include <sstream>
#include <utility>  // For pair

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// ExponentialFilter implementation
//==============================================================================
void
ExponentialFilter::from_xml(pugi::xml_node node)
{
    order_ = std::stoi(get_node_value(node, "order"));
    x_ = std::stod(get_node_value(node, "x"));
    y_ = std::stod(get_node_value(node, "y"));
    r_ = std::stod(get_node_value(node, "r"));
    exponent_ = std::stod(get_node_value(node, "exponent"));

    n_bins_ = order_ + 1;
}

void
ExponentialFilter::get_all_bins(const Particle* p, int estimator,
                            FilterMatch& match) const
{
  // Determine the normalized (r,theta) coordinates.
  double x = p->r().x - x_;
  double y = p->r().y - y_;
  double r = std::sqrt(x*x + y*y) / r_;

  if (r <= 1.0) {
    // Compute and return the Exponential weights.
    double exp[n_bins_];
    calc_exp(order_, exponent_, r, exp);
    for (int i = 0; i < n_bins_; i++) {
      match.bins_.push_back(i);
      match.weights_.push_back(exp[i]);
    }
  }
}

void
ExponentialFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);
  write_dataset(filter_group, "exponent", exponent_);
  write_dataset(filter_group, "x", x_);
  write_dataset(filter_group, "y", y_);
  write_dataset(filter_group, "r", r_);
}

std::string
ExponentialFilter::text_label(int bin) const
{
  std::stringstream out;
  out << "Exponential expansion, EXP" << std::to_string(bin);
  return out.str();
}

//==============================================================================
// C-API functions
//==============================================================================

std::pair<int, ExponentialFilter*>
check_exponential_filter(int32_t index)
{
  // Make sure this is a valid index to an allocated filter.
  int err = verify_filter(index);
  if (err) {
    return {err, nullptr};
  }

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<ExponentialFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Not an Exponential filter.");
    err = OPENMC_E_INVALID_TYPE;
  }
  return {err, filt};
}

extern "C" int
openmc_exponential_filter_get_order(int32_t index, int* order)
{
  // Check the filter.
  auto check_result = check_exponential_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the order.
  *order = filt->order_;
  return 0;
}

extern "C" int
openmc_exponential_filter_get_exponent(int32_t index, double* exponent)
{
  // Check the filter.
  auto check_result = check_exponential_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the exponent.
  *exponent = filt->exponent_;
  return 0;
}

extern "C" int
openmc_exponential_filter_get_params(int32_t index, double* x, double* y,
                                     double* r)
{
  // Check the filter.
  auto check_result = check_exponential_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Output the params.
  *x = filt->x_;
  *y = filt->y_;
  *r = filt->r_;
  return 0;
}

extern "C" int
openmc_exponential_filter_set_order(int32_t index, int order)
{
  // Check the filter.
  auto check_result = check_exponential_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  filt->order_ = order;
  filt->n_bins_ = order + 1;
  return 0;
}

extern "C" int
openmc_exponential_filter_set_exponent(int32_t index, double exponent)
{
  // Check the filter.
  auto check_result = check_exponential_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  filt->exponent_ = exponent;
  return 0;
}

extern "C" int
openmc_exponential_filter_set_params(int32_t index, const double* x,
                                 const double* y, const double* r)
{
  // Check the filter.
  auto check_result = check_exponential_filter(index);
  auto err = check_result.first;
  auto filt = check_result.second;
  if (err) return err;

  // Update the filter.
  if (x) filt->x_ = *x;
  if (y) filt->y_ = *y;
  if (r) filt->r_ = *r;
  return 0;
}

} // namespace openmc

