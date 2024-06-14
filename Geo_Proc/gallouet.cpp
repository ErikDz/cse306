#include "gallouet.h"
#include "fluid.h"
#include "rendering.h"
#include "liblbfgs/lbfgs.h" 

// Props to Milos Oundjian for big help on this again
// Debugging inspiration: https://github.com/ekaterina-borisova/graphics_cse306
void update_position(Vector &pos, Vector &vel, double dt) {
    pos = pos + dt * vel;
}

void apply_boundary_conditions(Vector &pos) {
    if (pos[0] < 0) pos[0] = -pos[0];
    if (pos[1] < 0) pos[1] = -pos[1];
    if (pos[0] >= 1) pos[0] = 2 - pos[0];
    if (pos[1] >= 1) pos[1] = 2 - pos[1];
}

void update_velocity(Vector &vel, const Vector &net_f, double dt, double m) {
    vel = vel + dt / m * net_f;
}

Vector calculate_spring_force(const Vector &centroid, const Vector &pos, double eps2) {
    return (centroid - pos) / eps2;
}

void apply_forces_and_update(std::vector<Vector> &pos_list, std::vector<Vector> &vel_list, std::vector<Polygon> &vor_cells, double m, double eps2, double dt) {
    Vector grav(0, -9.81);
    int np = pos_list.size();

    for (int i = 0; i < np; i++) {
        Vector spr_f = calculate_spring_force(vor_cells[i].calculate_centroid(), pos_list[i], eps2);
        Vector net_f = spr_f + m * grav;
        update_velocity(vel_list[i], net_f, dt, m);
        update_position(pos_list[i], vel_list[i], dt);
        apply_boundary_conditions(pos_list[i]);
    }
}

void gal_step(std::vector<Vector> &pos_list, std::vector<Vector> &vel_list, std::vector<double> &wt_list, int ts) {
    const double m = 200.0;
    const double eps = 0.004;
    const double dt = 0.02;
    double eps2 = std::pow(eps, 2);

    int np = pos_list.size();
    double obj_val;

    int opt_res = lbfgs(np + 1, &wt_list[0], &obj_val, eval_f, NULL, &pos_list[0], NULL);
    std::vector<Polygon> vor_cells = construct_fluid(&pos_list[0], &wt_list[0], np + 1);
    save_frame(vor_cells, "frames/frame_", ts);

    apply_forces_and_update(pos_list, vel_list, vor_cells, m, eps2, dt);
}