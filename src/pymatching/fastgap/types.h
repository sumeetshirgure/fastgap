// Copyright 2026 fastgap contributors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PYMATCHING2_FASTGAP_TYPES_H
#define PYMATCHING2_FASTGAP_TYPES_H

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace pm {
namespace fastgap {

/// All radii, distances and reduced costs are integers in PyMatching's internal units
/// (`w_int = 2 * round(w * c)`). See CLAUDE.md §2.1.
typedef int64_t dist_int;

/// "Unreachable" distance. Chosen so that sums of a few INF values cannot overflow int64.
constexpr dist_int DIST_INF = std::numeric_limits<int64_t>::max() / 8;

/// Mate codes for `DualState::mate`: a non-negative value is a defect index.
constexpr int32_t MATE_L = -1;
constexpr int32_t MATE_R = -2;

/// Boundary side codes stored per detector in `InteriorGraph::boundary_side`.
constexpr uint8_t SIDE_NONE = 0;
constexpr uint8_t SIDE_L = 1;  // boundary edge with obs' = 1
constexpr uint8_t SIDE_R = 2;  // boundary edge with obs' = 0

inline bool is_inf(dist_int d) {
    return d >= DIST_INF / 2;
}

/// Thrown when an invariant from CLAUDE.md §5.1 / §9 is violated. A violation is a bug.
struct InvariantViolation : std::logic_error {
    explicit InvariantViolation(const std::string& what) : std::logic_error("fastgap invariant violated: " + what) {
    }
};

#ifdef FASTGAP_CHECK
constexpr bool CHECKS_DEFAULT_ON = true;
#else
constexpr bool CHECKS_DEFAULT_ON = false;
#endif

}  // namespace fastgap
}  // namespace pm

#endif  // PYMATCHING2_FASTGAP_TYPES_H
