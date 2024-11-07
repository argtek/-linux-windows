/*
 * KF-GINS: An EKF-Based GNSS/INS Integrated Navigation System
 *
 * Copyright (C) 2022 i2Nav Group, Wuhan University
 *
 *     Author : Liqiang Wang
 *    Contact : wlq@whu.edu.cn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef GI_ENGINE_H
#define GI_ENGINE_H

#define IN
#define OUT


#include <Eigen/Dense>
#include <vector>

#include "common/types.h"

#include "kf_gins_types.h"

class GIEngine
{

public:
    explicit GIEngine(GINSOptions& options);

    ~GIEngine() = default;

    /**
     * @brief 添加新的IMU数据，(不)补偿IMU误差
     *        add new imudata, do (not) compensate imu error
     * @param [in] imu        新的IMU原始数据
     *                        new raw imudata
     * @param [in] compensate 是否补偿IMU误差
     *                        if compensate imu error to new imudata
     * */
    void addImuData(const IMU& imu, bool compensate = false)
    {

        m_imuPre_ = m_imuCur_;
        m_imuCur_ = imu;

        if (compensate)
        {
            imuCompensate(m_imuCur_);
        }
    }

    /**
     * @brief 添加新的GNSS数据
     *        add new gnssdata
     * @param [in] gnss 新的GNSS数据
     *                  new gnssdata
     * */
    void addGnssData(const GNSS& gnss)
    {

        m_gnssData_ = gnss;
        // 暂不进行数据有效性检查，GNSS数据默认有效
        // do not check the validity of gnssdata, the gnssdata is valid by default
        m_gnssData_.isvalid = true;
    }

    /**
     * @brief 处理新的IMU数据
     *        process new imudata
     * */
    void newImuDataProcess();

    /**
     * @brief 内插增量形式的IMU数据到指定时刻
     *        interpolate incremental imudata to given timestamp
     * @param [in]     imu1      前一时刻IMU数据
     *                           the previous imudata
     * @param [in,out] imu2      当前时刻IMU数据
     *                           the current imudata
     * @param [in]     timestamp 给定内插到的时刻
     *                           given interpolate timestamp
     * @param [in,out] midimu    输出内插时刻的IMU数据
     *                           output imudata at given timestamp
     * */
    static void imuInterpolate(IN const IMU& imu1, IN OUT IMU& imu2, IN const double timestamp, IN OUT IMU& midimu)
    {

        if (imu1.time > timestamp || imu2.time < timestamp)
        {
            return;
        }

        double lamda = (timestamp - imu1.time) / (imu2.time - imu1.time);

        midimu.time = timestamp;
        midimu.dtheta = imu2.dtheta * lamda;
        midimu.dvel = imu2.dvel * lamda;
        midimu.dt = timestamp - imu1.time;

        imu2.dtheta = imu2.dtheta - midimu.dtheta;
        imu2.dvel = imu2.dvel - midimu.dvel;
        imu2.dt = imu2.dt - midimu.dt;
    }

    /**
     * @brief 获取当前时间
     *        get current time
     * */
    double timestamp() const
    {
        return m_timeStamp_;
    }

    /**
     * @brief 获取当前IMU状态
     *        get current navigation state
     * */
    NavState getNavState();

    /**
     * @brief 获取当前状态协方差
     *        get current state covariance
     * */
    Eigen::MatrixXd getCovariance()
    {
        return m_Cov_;
    }

private:
    /**
     * @brief 初始化系统状态和协方差
     *        initialize state and state covariance
     * @param [in] initstate     初始状态
     *                           initial state
     * @param [in] initstate_std 初始状态标准差
     *                           initial state std
     * */
    void initialize(IN const NavState& initstate, IN const NavState& initstate_std);

    /**
     * @brief 当前IMU误差补偿到IMU数据中
     *        componsate imu error to the imudata
     * @param [in,out] imu 需要补偿的IMU数据
     *                     imudata to be compensated
     * */
    void imuCompensate(IN OUT IMU& imu);

    /**
     * @brief 判断是否需要更新,以及更新哪一时刻系统状态
     *        determine if we should do upate and which navstate to update
     * @param [in] imutime1   上一IMU状态时间
     *                        the last state time
     * @param [in] imutime2   当前IMU状态时间
     *                        the current state time
     * @param [in] updatetime 状态更新的时间
     *                        time to update state
     * @return 0: 不需要更新
     *            donot need update
     *         1: 需要更新上一IMU状态
     *            update the last navstate
     *         2: 需要更新当前IMU状态
     *            update the current navstate
     *         3: 需要将IMU进行内插到状态更新时间
     *            need interpolate imudata to updatetime
     * */
    int isToUpdate(double imutime1, double imutime2, double updatetime) const;

    /**
     * @brief 进行INS状态更新(IMU机械编排算法), 并计算IMU状态转移矩阵和噪声阵
     *        do INS state update(INS mechanization), and compute state transition matrix and noise matrix
     * @param [in,out] imupre 前一时刻IMU数据
     *                        imudata at the previous epoch
     * @param [in,out] imucur 当前时刻IMU数据
     *                        imudata at the current epoch
     * */
    void insPropagation(IN OUT IMU& imupre, IN OUT IMU& imucur);

    /**
     * @brief 使用GNSS位置观测更新系统状态
     *        update state using gnss position
     * @param [in,out] gnssdata
     * */
    void gnssUpdate(IN OUT GNSS& gnssdata);

    /**
     * @brief Kalman 预测,
     *        Kalman Filter Predict process
     * @param [in,out] Phi 状态转移矩阵
     *                     state transition matrix
     * @param [in,out] Qd  传播噪声矩阵
     *                     propagation noise matrix
     * */
    void EKFPredict(IN OUT Eigen::MatrixXd& Phi, IN OUT Eigen::MatrixXd& Qd);

    /**
     * @brief Kalman 更新
     *        Kalman Filter Update process
     * @param [in] dz 观测新息
     *                measurement innovation
     * @param [in] H  观测矩阵
     *                measurement matrix
     * @param [in] R  观测噪声阵
     *                measurement noise matrix
     * */
    void EKFUpdate(IN Eigen::MatrixXd& dz, IN Eigen::MatrixXd& H, IN Eigen::MatrixXd& R);

    /**
     * @brief 反馈误差状态到当前状态
     *        feedback error state to the current state
     * */
    void stateFeedback();

    /**
     * @brief 检查协方差对角线元素是否都为正
     *        Check if covariance diagonal elements are all positive
     * */
    void checkCov()
    {

        for (int i = 0; i < m_constRANK; i++)
        {
            if (m_Cov_(i, i) < 0)
            {
                std::cout << "Covariance is negative at " << std::setprecision(10) << m_timeStamp_ << " !" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }

private:
    GINSOptions m_options_;

    double m_timeStamp_;

    // 更新时间对齐误差，IMU状态和观测信息误差小于它则认为两者对齐
    // updata time align error
    const double m_constTIME_ALIGN_ERR = 0.001;

    // IMU和GNSS原始数据
    // raw imudata and gnssdata
    IMU m_imuPre_;
    IMU m_imuCur_;
    GNSS m_gnssData_;

    // IMU状态（位置、速度、姿态和IMU误差）
    // imu state (position, velocity, attitude and imu error)
    PVA m_pvaCur_;
    PVA m_pvaPre_;
    ImuError m_imuError_;

    // Kalman滤波相关
    // ekf variables
    Eigen::MatrixXd m_Cov_;
    Eigen::MatrixXd m_Qc_;
    Eigen::MatrixXd m_dx_;

    const int m_constRANK = 21;
    const int m_constNOISERANK = 18;

    // 状态ID和噪声ID
    // state ID and noise ID
    enum StateID 
    { 
        // 状态ID
        P_ID = 0,  // 位置状态
        V_ID = 3, 
        PHI_ID = 6, // 状态转移矩阵通常用希腊字母 Φ（phi）表示
        BG_ID = 9,  // 零偏
        BA_ID = 12, 
        SG_ID = 15, // 比例因子
        SA_ID = 18 
    };

    enum NoiseID 
    { 
        // 噪声ID
        VRW_ID = 0, // 随机游走噪声
        ARW_ID = 3, 
        BGSTD_ID = 6,  // 零偏标准差
        BASTD_ID = 9, 
        SGSTD_ID = 12,  // 比例因子标准差
        SASTD_ID = 15 
    };
};

#endif // GI_ENGINE_H
