#include <windows.h>
#include <vector>
#include <C:\Users\leesp\Downloads\boost_1_88_0\boost\numeric\odeint.hpp>
#include <limits>

#ifdef __cplusplus
extern "C" {
#endif

#define DLL_EXPORT __declspec(dllexport)

struct Vector3 {
    double x, y, z;
    Vector3(double _x = 0.0, double _y = 0.0, double _z = 0.0) : x(_x), y(_y), z(_z) {}
};

Vector3 cross(const Vector3& a, const Vector3& b) {
    return Vector3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

Vector3 inertia_tensor_multiply(const Vector3& J_diag, const Vector3& v) {
    return Vector3(J_diag.x * v.x, J_diag.y * v.y, J_diag.z * v.z);
}

struct RigidBodySystem {
    Vector3 J_diag;
    Vector3 r;
    double m;
    Vector3 g;

    RigidBodySystem(const Vector3& _J_diag, const Vector3& _r, double _m, const Vector3& _g)
        : J_diag(_J_diag), r(_r), m(_m), g(_g) {}

    void operator()(const std::vector<double>& state, std::vector<double>& dstate_dt, double /* t */) const {
        Vector3 w(state[0], state[1], state[2]);
        Vector3 force(m * g.x, m * g.y, m * g.z);
        Vector3 torque = cross(r, force);
        Vector3 J_omega = inertia_tensor_multiply(J_diag, w);
        Vector3 cross_term = cross(w, J_omega);
        Vector3 J_dot_omega = Vector3(
            torque.x - cross_term.x,
            torque.y - cross_term.y,
            torque.z - cross_term.z
        );

        dstate_dt[0] = J_dot_omega.x / J_diag.x;
        dstate_dt[1] = J_dot_omega.y / J_diag.y;
        dstate_dt[2] = J_dot_omega.z / J_diag.z;

        double R[3][3] = {
            {state[3], state[4], state[5]},
            {state[6], state[7], state[8]},
            {state[9], state[10], state[11]}
        };

        double skew_omega[3][3] = {
            {   0.0, -w.z,  w.y},
            {  w.z,   0.0, -w.x},
            { -w.y,  w.x,   0.0}
        };

        double dR_dt[3][3] = { {0.0} };
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    dR_dt[i][j] += R[i][k] * skew_omega[k][j];
                }
            }
        }

        dstate_dt[3] = dR_dt[0][0];
        dstate_dt[4] = dR_dt[0][1];
        dstate_dt[5] = dR_dt[0][2];
        dstate_dt[6] = dR_dt[1][0];
        dstate_dt[7] = dR_dt[1][1];
        dstate_dt[8] = dR_dt[1][2];
        dstate_dt[9] = dR_dt[2][0];
        dstate_dt[10] = dR_dt[2][1];
        dstate_dt[11] = dR_dt[2][2];
    }
};

DLL_EXPORT double SolveRigidBodyODE(
    double J_xx, double J_yy, double J_zz,
    double m,
    double r_x, double r_y, double r_z,
    double w0_x, double w0_y, double w0_z,
    double R0_11, double R0_12, double R0_13,
    double R0_21, double R0_22, double R0_23,
    double R0_31, double R0_32, double R0_33,
    double g1, double g2, double g3,
    double* w_array,
    double* R_array
) {
    if (!w_array || !R_array ||
        J_xx <= 0 || J_yy <= 0 || J_zz <= 0 || m <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double dt = 0.001;
    std::vector<double> state = {
        w0_x, w0_y, w0_z,
        R0_11, R0_12, R0_13,
        R0_21, R0_22, R0_23,
        R0_31, R0_32, R0_33
    };
    std::vector<double> solution(12);

    Vector3 J_diag(J_xx, J_yy, J_zz);
    Vector3 r(r_x, r_y, r_z);
    Vector3 g(g1, g2, g3);
    RigidBodySystem system(J_diag, r, m, g);

    try {
        boost::numeric::odeint::integrate_const(
            boost::numeric::odeint::runge_kutta4<std::vector<double>>(),
            system,
            state,
            0.0, // t0 寫死為 0
            dt,
            dt,
            [&](const std::vector<double>& state_val, double t_val) {
                solution = state_val;
            }
        );

        w_array[0] = solution[0]; // wx
        w_array[1] = solution[1]; // wy
        w_array[2] = solution[2]; // wz
        for (int j = 0; j < 9; ++j) {
            R_array[j] = solution[3 + j];
        }
        return solution[0]; // 返回最終 wx
    }
    catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

BOOL APIENTRY DllMain(HMODULE hModule, DWORD ul_reason_for_call, LPVOID lpReserved) {
    switch (ul_reason_for_call) {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

#ifdef __cplusplus
}
#endif