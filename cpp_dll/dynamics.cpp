#include <windows.h>
#include <vector>
#include <C:\Users\leesp\Downloads\boost_1_88_0\boost\numeric\odeint.hpp>

#ifdef __cplusplus
extern "C" {
#endif

// 導出函數的宏
#define DLL_EXPORT __declspec(dllexport)

// 3D 向量結構
struct Vector3 {
    double x, y, z;
    Vector3(double _x = 0.0, double _y = 0.0, double _z = 0.0) : x(_x), y(_y), z(_z) {}
};

// 叉積：a × b
Vector3 cross(const Vector3& a, const Vector3& b) {
    return Vector3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

// 慣性張量與向量的乘法：J * v (假設 J 為對角矩陣)
Vector3 inertia_tensor_multiply(const Vector3& J_diag, const Vector3& v) {
    return Vector3(J_diag.x * v.x, J_diag.y * v.y, J_diag.z * v.z);
}

// 微分方程系統：角速度和姿態
struct RigidBodySystem {
    Vector3 J_diag; // 慣性張量對角元素
    Vector3 r;      // 重力作用點
    double m;       // 質量
    Vector3 g;      // 重力加速度

    RigidBodySystem(const Vector3& _J_diag, const Vector3& _r, double _m, const Vector3& _g)
        : J_diag(_J_diag), r(_r), m(_m), g(_g) {}

    void operator()(const std::vector<double>& state, std::vector<double>& dstate_dt, double /* t */) const {
        // 狀態向量：state = [wx, wy, wz, R11, R12, R13, R21, R22, R23, R31, R32, R33]
        Vector3 w(state[0], state[1], state[2]);

        // 1. 角速度動態：J * dot(omega) + omega × (J * omega) = r × mg
        Vector3 torque = cross(r, Vector3(0.0, 0.0, -m * 9.81));
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

        // 2. 姿態動態：dR/dt = R * skew(omega)
        // 提取 R 矩陣
        double R[3][3] = {
            {state[3], state[4], state[5]},
            {state[6], state[7], state[8]},
            {state[9], state[10], state[11]}
        };

        // 角速度的斜對稱矩陣
        double skew_omega[3][3] = {
            {   0.0, -w.z,  w.y},
            {  w.z,   0.0, -w.x},
            { -w.y,  w.x,   0.0}
        };

        // 計算 dR/dt = R * skew(omega)
        double dR_dt[3][3] = { {0.0} };
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    dR_dt[i][j] += R[i][k] * skew_omega[k][j];
                }
            }
        }

        // 填充 dstate_dt 的姿態部分
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

// 求解 ODE 並將結果寫入輸出陣列
DLL_EXPORT int32_t SolveRigidBodyODE(
    double J_xx, double J_yy, double J_zz, // 慣性張量對角元素
    double m,                              // 質量
    double r_x, double r_y, double r_z,    // 重力作用點
    double w0_x, double w0_y, double w0_z, // 初始角速度
    double R0_11, double R0_12, double R0_13, // 初始姿態矩陣
    double R0_21, double R0_22, double R0_23,
    double R0_31, double R0_32, double R0_33,
    double t0, double t_end, int32_t n,    // 時間範圍和步數
    double* time_array,                    // 輸出：時間陣列
    double* wx_array, double* wy_array, double* wz_array, // 輸出：角速度
    double* R_array // 輸出：姿態矩陣（n * 9 元素）
) {
    // 參數驗證
    if (n <= 0 || !time_array || !wx_array || !wy_array || !wz_array || !R_array ||
        t_end <= t0 || J_xx <= 0 || J_yy <= 0 || J_zz <= 0 || m <= 0) {
        return -1; // 錯誤碼：無效參數
    }

    // 設置時間點
    double dt = (t_end - t0) / (n - 1);
    std::vector<double> t(n);
    for (int32_t i = 0; i < n; ++i) {
        t[i] = t0 + i * dt;
    }

    // 初始條件：omega 和 R
    std::vector<double> state = {
        w0_x, w0_y, w0_z, // omega
        R0_11, R0_12, R0_13, // R 矩陣（按行展開）
        R0_21, R0_22, R0_23,
        R0_31, R0_32, R0_33
    };
    std::vector<std::vector<double>> solution(n, std::vector<double>(12));

    // 設置系統參數
    Vector3 J_diag(J_xx, J_yy, J_zz);
    Vector3 r(r_x, r_y, r_z);
    Vector3 g(0.0, 0.0, -9.81);
    RigidBodySystem system(J_diag, r, m, g);

    // 使用 odeint 求解
    try {
        boost::numeric::odeint::integrate_const(
            boost::numeric::odeint::runge_kutta4<std::vector<double>>(),
            system,
            state,
            t0,
            t_end,
            dt,
            [&](const std::vector<double>& state_val, double t_val) {
                int32_t idx = static_cast<int32_t>((t_val - t0) / dt + 0.5);
                if (idx >= 0 && idx < n) {
                    solution[idx] = state_val;
                }
            }
        );

        // 填充輸出陣列
        for (int32_t i = 0; i < n; ++i) {
            time_array[i] = t[i];
            wx_array[i] = solution[i][0];
            wy_array[i] = solution[i][1];
            wz_array[i] = solution[i][2];
            // R 矩陣：9 個元素
            for (int j = 0; j < 9; ++j) {
                R_array[i * 9 + j] = solution[i][3 + j];
            }
        }
        return 0; // 成功
    }
    catch (...) {
        return -2; // 錯誤碼：求解失敗
    }
}

// DLL 初始化
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