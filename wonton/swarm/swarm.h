/*
 * This file is part of the Ristra wonton project.
 * Please see the license file at the root of this repository, or at:
 *  https://github.com/laristra/wonton/blob/master/LICENSE
 */
#ifndef WONTON_SWARM_H_
#define WONTON_SWARM_H_

#include <vector>
#include <memory>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <iostream>
#include <random>

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

namespace Wonton {

/**
 * @class Swarm
 *
 * @brief An effective "mesh" class for a collection disconnected points (particles).
 *
 * @tparam dim: the dimension of the problem.
 */
template<int dim>
class Swarm {
public:

  /**
   * @brief Create an empty swarm.
   *
   */
  Swarm() = default;

  /**
   * @brief Create a particle field from a given list.
   *
   * @param points: the given list of points.
   */
  explicit Swarm(Wonton::vector<Wonton::Point<dim>> const& points)
    : points_(points),
      num_local_points_(points_.size())
  {}

#ifdef WONTON_ENABLE_THRUST
  /**
   * @brief Create a particle field from a given list.
   *
   * @param points: the given list of points.
   */
  explicit Swarm(std::vector<Wonton::Point<dim>> const& points)
    : points_(points),
      num_local_points_(points_.size())
  {}

  /**
   * @brief Initialize from a given list of points.
   *
   * @param points: the given list of points.
   */
  inline void init(std::vector<Wonton::Point<dim>> const& points) {
    num_local_points_ = points.size();
    points_.resize(num_local_points_);
    std::copy(points.begin(), points.end(), points_.begin());
  }
#endif

  /**
   * @brief Initialize from a given list of points.
   *
   * @param points: the given list of points.
   */
  inline void init(Wonton::vector<Wonton::Point<dim>> const& points) {
    num_local_points_ = points.size();
    points_.resize(num_local_points_);
    std::copy(points.begin(), points.end(), points_.begin());
  }

  /**
   * @brief Resize internal list of points.
   *
   * @param new_size: the new size
   */
  inline void resize(int new_size) {
    assert(new_size > 0);
    num_local_points_ = new_size;
    points_.resize(new_size);
  }

  /**
   * @brief Assign the given point at the given location.
   *
   * @param id: index of the point.
   * @param p: coordinates of the point.
   */
  inline void assign(int id, Wonton::Point<dim> const& p) {
    assert(id >= 0 and id < num_local_points_);
    points_[id] = p;
  }

  /**
   * @brief Generate a random or uniform particle field.
   *
   * @param num_particles: number of points.
   * @param p_min: lower bounds on points coordinates.
   * @param p_max: upper bounds on points coordinates.
   * @param distribution: the kind of random numbers distribution.
   * @param seed: a seed used by the random number generator.
   */
  Swarm(int num_particles, int distribution,
        unsigned user_seed = 0,
        double x_min = 0.0, double x_max = 0.0,
        double y_min = 0.0, double y_max = 0.0,
        double z_min = 0.0, double z_max = 0.0);

  /**
   * @brief Create a particle field from an arbitary mesh wrapper.
   *
   * @tparam Mesh: type of mesh wrapper.
   * @param mesh: a mesh.
   * @param entity: the entity kind.
   */
  template<typename Mesh>
  Swarm(Mesh const& mesh, Wonton::Entity_kind entity);

  /**
   * @brief Create a particle field from a set of mesh wrappers.
   *
   * @tparam Mesh: type of mesh wrappers.
   * @param meshes: a set of meshes.
   * @param entity: the entity kind.
   */
  template<typename Mesh>
  Swarm(std::vector<Mesh*> const& meshes, Wonton::Entity_kind entity);

  /**
   * @brief Get the dimension.
   *
   * @return the dimension.
   */
  int space_dimension() const { return dim; }

  /**
   * @brief Get owned particles count.
   *
   * @return number of owned particles.
   */
  int num_owned_particles() const { return num_local_points_; }

  /**
   * @brief Get ghost particles count.
   *
   * @return number of ghost particles.
   */
  int num_ghost_particles() const { return points_.size() - num_local_points_; }

  /**
   * @brief Get particles count
   *
   * @param type: entity type (owned, ghost, all).
   * @return the particle count for the given entity type.
   */
  int num_particles(Wonton::Entity_type type = Wonton::ALL) const {
    switch (type) {
      case Wonton::PARALLEL_OWNED: return num_owned_particles();
      case Wonton::PARALLEL_GHOST: return num_ghost_particles();
      default: return num_owned_particles() + num_ghost_particles();
    }
  }

  /**
   * @brief Get the coordinates of the particle.
   *
   * @param index: index of particle.
   * @return its coordinates.
   */
  Wonton::Point<dim> get_particle_coordinates(int index) const {
    Wonton::Point<dim> point = points_[index];
    return point;
  }

  /**
   * @brief Get an iterator at first element of the particle field.
   *
   * @param kind: entity kind (kept for consistency with mesh iterators)
   * @param type: entity type (kept for consistency with mesh iterators)
   * @return an iterator pointing to the first element of the particle field.
   */
  counting_iterator begin(Entity_kind const kind = Wonton::PARTICLE,
                          Entity_type const type = Wonton::ALL) const {

    assert(kind == Wonton::PARTICLE);
    return make_counting_iterator(0);
  }

  /**
   * @brief Get an iterator at last element of the particle field.
   *
   * @param kind: entity kind (kept for consistency with mesh iterators)
   * @param type: entity type (kept for consistency with mesh iterators)
   * @return an iterator pointing to the last element of the particle field.
   */
  counting_iterator end(Entity_kind const kind = Wonton::PARTICLE,
                        Entity_type const type = Wonton::ALL) const {

    assert(kind == Wonton::PARTICLE);
    return (make_counting_iterator(0) + num_particles(type));
  }

  /**
   * @brief Extend the particle field.
   *
   */
  void extend_particle_list(std::vector<Wonton::Point<dim>> const& new_pts) {
    points_.insert(points_.end(), new_pts.begin(), new_pts.end());
  }

private:
  /** the centers of the particles */
  Wonton::vector<Wonton::Point<dim>> points_;

  /** the number of owned particles in the swarm */
  int num_local_points_ = 0;
};

/* -------------------------------------------------------------------------- */
template<>
inline Swarm<1>::Swarm(int num_particles, int distribution, unsigned user_seed,
                       double x_min, double x_max,
                       double /* unused */, double /* unused */,
                       double /* unused */, double /* unused */) {

  assert(num_particles > 0);
  assert(0 <= distribution and distribution <= 2);

  // resize field
  points_.resize(num_particles);
  num_local_points_ = num_particles;

  // set the random engine and generator
  std::random_device device;
  std::mt19937 engine { user_seed ? user_seed : device() };
  std::uniform_real_distribution<double> generator(0.0, 1.0);

  // update coordinates
  if (distribution == 0) {
    for (auto&& current : points_)
      current = Wonton::Point<1>(x_min + (x_max - x_min) * generator(engine));

  } else {
    double const h = (x_max - x_min) / (num_particles - 1);
    for (int i = 0; i < num_particles; i++)
      points_[i] = Wonton::Point<1>(x_min + i * h);

    if (distribution == 2) {
      for (int i = 0; i < num_particles; i++) {
        Wonton::Point<1> current = points_[i];
        current[0] += 0.25 * h * (2 * generator(engine) - 1);
        points_[i] = Wonton::Point<1>(std::max(x_min, std::min(x_max, current[0])));
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline Swarm<2>::Swarm(int num_particles, int distribution, unsigned user_seed,
                       double x_min, double x_max,
                       double y_min, double y_max,
                       double /* unused */, double /* unused */) {

  assert(num_particles > 0);
  assert(0 <= distribution and distribution <= 2);

  // set the random engine and generator
  std::random_device device;
  std::mt19937 engine { user_seed ? user_seed : device() };
  std::uniform_real_distribution<double> generator(0.0, 1.0);

  if (distribution == 0) {
    // resize field and update coordinates
    num_local_points_ = num_particles;
    points_.resize(num_local_points_);

    for (auto&& current : points_) {
      // point coordinates are not always initialized in order
      // so enforce random number picking sequence
      double const noise[] = { generator(engine),
                               generator(engine) };

      current = Wonton::Point<2>(x_min + (x_max - x_min) * noise[0],
                                 y_min + (y_max - y_min) * noise[1]);
    }
  } else {
    // resize field
    int const num_per_dim = std::floor(std::sqrt(num_particles));
    num_local_points_ = num_per_dim * num_per_dim;
    points_.resize(num_local_points_);
    double const hx = (x_max - x_min) / (num_per_dim - 1);
    double const hy = (y_max - y_min) / (num_per_dim - 1);

    // update coordinates
    int k = 0;
    for (int i = 0; i < num_per_dim; i++)
      for (int j = 0; j < num_per_dim; ++j)
        points_[k++] = Wonton::Point<2>(x_min + i * hx, y_min + j * hy);

    if (distribution == 2) {
      for (auto&& current : points_) {
        // point coordinates are not always initialized in order
        // so enforce random number picking sequence
        double const noise[] = { generator(engine),
                                 generator(engine) };

        Wonton::Point<2> copy = current;
        copy[0] += 0.25 * hx * (2 * noise[0] - 1);
        copy[1] += 0.25 * hy * (2 * noise[1] - 1);
        current = Wonton::Point<2>(std::max(x_min, std::min(x_max, copy[0])),
                                   std::max(y_min, std::min(y_max, copy[1])));
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline Swarm<3>::Swarm(int num_particles, int distribution, unsigned user_seed,
                       double x_min, double x_max,
                       double y_min, double y_max,
                       double z_min, double z_max) {

  assert(num_particles > 0);
  assert(0 <= distribution and distribution <= 2);

  // set the random engine and generator
  std::random_device device;
  std::mt19937 engine { user_seed ? user_seed : device() };
  std::uniform_real_distribution<double> generator(0.0, 1.0);

  if (distribution == 0) {
    // resize field and update coordinates
    num_local_points_ = num_particles;
    points_.resize(num_local_points_);

    for (int i = 0; i < num_particles; i++) {
      // point coordinates are not always initialized in order
      // so enforce random number picking sequence
      double const noise[] = { generator(engine),
                               generator(engine),
                               generator(engine) };

      points_[i] = Wonton::Point<3>(x_min + (x_max - x_min) * noise[0],
                                    y_min + (y_max - y_min) * noise[1],
                                    z_min + (z_max - z_min) * noise[2]);
    }

  } else {
    auto const cubic_root = std::pow(num_particles, 1./3.);
    auto const num_per_dim = static_cast<int>(std::ceil(cubic_root));
    auto const hx = (x_max - x_min) / (cubic_root - 1);
    auto const hy = (y_max - y_min) / (cubic_root - 1);
    auto const hz = (z_max - z_min) / (cubic_root - 1);

    // resize field
    num_local_points_ = num_per_dim * num_per_dim * num_per_dim;
    points_.resize(num_local_points_);

    // update coordinates
    int n = 0;
    for (int i = 0; i < num_per_dim; i++)
      for (int j = 0; j < num_per_dim; ++j)
        for (int k = 0; k < num_per_dim; ++k)
          points_[n++] = Wonton::Point<3>(x_min + i * hx,
                                          y_min + j * hy,
                                          z_min + k * hz);

    if (distribution == 2) {
      for (auto&& current : points_) {
        // point coordinates are not always initialized in order
        // so enforce random number picking sequence
        double const noise[] = { generator(engine),
                                 generator(engine),
                                 generator(engine) };

        Wonton::Point<3> copy = current;
        copy[0] += 0.25 * hx * (2 * noise[0] - 1);
        copy[1] += 0.25 * hy * (2 * noise[1] - 1);
        copy[2] += 0.25 * hz * (2 * noise[2] - 1);
        current = Wonton::Point<3>(std::max(x_min, std::min(x_max, copy[0])),
                                   std::max(y_min, std::min(y_max, copy[1])),
                                   std::max(z_min, std::min(z_max, copy[2])));
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<int dim>
template<typename Mesh>
Swarm<dim>::Swarm(Mesh const& mesh, Wonton::Entity_kind entity) {
  switch (entity) {
    case Wonton::NODE: {
      num_local_points_ = mesh.num_owned_nodes();
      points_.resize(num_local_points_);
      Wonton::Point<dim> coord;
      for (int i = 0; i < num_local_points_; ++i) {
        mesh.node_get_coordinates(i, &coord);
        points_[i] = coord;
      }
    } break;
    case Wonton::CELL: {
      num_local_points_ = mesh.num_owned_cells();
      points_.resize(num_local_points_);
      Wonton::Point<dim> centroid;
      for (int i = 0; i < num_local_points_; ++i) {
        mesh.cell_centroid(i, &centroid);
        points_[i] = centroid;
      }
    } break;
    default: throw std::runtime_error("unsupported entity");
  }
}

/* -------------------------------------------------------------------------- */
template<int dim>
template<typename Mesh>
Swarm<dim>::Swarm(std::vector<Mesh*> const& meshes, Wonton::Entity_kind entity) {

  // resize particle field
  int n = 0;
  switch (entity) {
    case Wonton::NODE: for (auto&& mesh : meshes) n += mesh->num_owned_nodes(); break;
    case Wonton::CELL: for (auto&& mesh : meshes) n += mesh->num_owned_cells(); break;
    default: throw std::runtime_error("unsupported entity");
  }

  num_local_points_ = n;
  points_.resize(num_local_points_);

  n = 0;
  switch (entity) {
    case Wonton::NODE: {
      for (auto&& mesh : meshes) {
        int const num_nodes = mesh->num_owned_nodes();
        Wonton::Point<dim> coord;
        for (int i = 0; i < num_nodes; i++) {
          mesh->node_get_coordinates(i, &coord);
          points_[n++] = coord;
        }
      }
    } break;
    case Wonton::CELL: {
      for (auto&& mesh : meshes) {
        int const num_cells = mesh->num_owned_cells();
        Wonton::Point<dim> centroid;
        for (int i = 0; i < num_cells; i++) {
          mesh->cell_centroid(i, &centroid);
          points_[n++] = centroid;
        }
      }
    } break;
    default: break;
  }
}

/* -------------------------------------------------------------------------- */
}  // namespace Wonton

#endif  // WONTON_SWARM_H_