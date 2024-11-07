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

#include "common/earth.h"
#include "common/rotation.h"

#include "gi_engine.h"
#include "insmech.h"

GIEngine::GIEngine(GINSOptions& options)
{

    this->m_options_ = options;
    m_options_.print_options();
    m_timeStamp_ = 0;

    // 设置协方差矩阵，系统噪声阵和系统误差状态矩阵大小
    // resize covariance matrix, system noise matrix, and system error state matrix
    m_Cov_.resize(m_constRANK, m_constRANK);
    m_Qc_.resize(m_constNOISERANK, m_constNOISERANK);
    m_dx_.resize(m_constRANK, 1);
    m_Cov_.setZero();
    m_Qc_.setZero();
    m_dx_.setZero();

    // 初始化系统噪声阵
    // initialize noise matrix
    auto imunoise = m_options_.imunoise;
    m_Qc_.block(ARW_ID, ARW_ID, 3, 3) = imunoise.gyr_arw.cwiseProduct(imunoise.gyr_arw).asDiagonal();
    m_Qc_.block(VRW_ID, VRW_ID, 3, 3) = imunoise.acc_vrw.cwiseProduct(imunoise.acc_vrw).asDiagonal();
    m_Qc_.block(BGSTD_ID, BGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrbias_std.cwiseProduct(imunoise.gyrbias_std).asDiagonal();
    m_Qc_.block(BASTD_ID, BASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accbias_std.cwiseProduct(imunoise.accbias_std).asDiagonal();
    m_Qc_.block(SGSTD_ID, SGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrscale_std.cwiseProduct(imunoise.gyrscale_std).asDiagonal();
    m_Qc_.block(SASTD_ID, SASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accscale_std.cwiseProduct(imunoise.accscale_std).asDiagonal();

    // 设置系统状态(位置、速度、姿态和IMU误差)初值和初始协方差
    // set initial state (position, velocity, attitude and IMU error) and covariance
    initialize(m_options_.initstate, m_options_.initstate_std);
}

void GIEngine::initialize(const NavState& initstate, const NavState& initstate_std)
{

    // 需要确定初始姿态，初始姿态是后续位置的基础
    {
        // 初始化位置、速度、姿态
        // initialize position, velocity and attitude
        m_pvaCur_.pos = initstate.pos;
        m_pvaCur_.vel = initstate.vel;
        m_pvaCur_.att.euler = initstate.euler;
        m_pvaCur_.att.cbn = Rotation::euler2matrix(m_pvaCur_.att.euler);
        m_pvaCur_.att.qbn = Rotation::euler2quaternion(m_pvaCur_.att.euler);
    }

    {
        // 噪声：零偏，比例因子，这个和芯片本身相关，有芯片厂商提供
        // 初始化IMU误差
        // initialize imu error
        m_imuError_ = initstate.imuerror;
    }

    // PVA：载体的姿态描述，就是最终要获得的姿态
    {
        // 给上一时刻状态赋同样的初值
        // set the same value to the previous state
        m_pvaPre_ = m_pvaCur_;
    }

    // 初始化协方差（协方差，两个样本x,y的相关程度，不一定是线性关系，可能是正负相关关系）
    // initialize covariance （co：协同，variance：变量）
    ImuError imuerror_std = initstate_std.imuerror; // 方差，标准差。对样本和均值之间差距的描述

    // block：提取矩阵里面子矩阵的函数
    /*
    * MatrixType block(int startRow, int startCol, int blockRows, int blockCols)
    * startRow：子矩阵的起始行索引。
    * startCol：子矩阵的起始列索引。
    * blockRows：子矩阵的行数。
    * blockCols：子矩阵的列数。
    */
    // cwiseProduct：对应位置矩阵元素相乘，输出到对应位置，不是矩阵乘法
    // asDiagonal 是一个用于创建对角矩阵的函数。通过该函数，可以将一个向量转换成对应的对角矩阵，其中向量的元素将成为对角矩阵的对角线元素，而其他位置的元素则为零
    // pos：向量，nitstate_std.pos.cwiseProduct(initstate_std.pos)：自身相乘，就是向量内每个元素的平方，将nitstate_std这个方差平方后转换成为方差
    {
        // test
        int CovRow=m_Cov_.rows();
        int CovCol=m_Cov_.cols();
    }
    // 21维协方差矩阵，描述了所有状态因素，每种状态都是用矩阵描述，向量也转化成为矩阵描述
    // 0，0，3，3：位置矩阵
    m_Cov_.block(P_ID, P_ID, 3, 3) = initstate_std.pos.cwiseProduct(initstate_std.pos).asDiagonal();

    // 3，3，6，6：速度矩阵
    m_Cov_.block(V_ID, V_ID, 3, 3) = initstate_std.vel.cwiseProduct(initstate_std.vel).asDiagonal();

    // 6，6，9，9：Phi:状态转移矩阵
    m_Cov_.block(PHI_ID, PHI_ID, 3, 3) = initstate_std.euler.cwiseProduct(initstate_std.euler).asDiagonal();

    // 9，9，12，12：陀螺仪零偏矩阵
    m_Cov_.block(BG_ID, BG_ID, 3, 3) = imuerror_std.gyrbias.cwiseProduct(imuerror_std.gyrbias).asDiagonal();

    // 12，12，15，15：加速度计零偏矩阵
    m_Cov_.block(BA_ID, BA_ID, 3, 3) = imuerror_std.accbias.cwiseProduct(imuerror_std.accbias).asDiagonal();

    // 15，15，18，18：陀螺仪比例因子
    m_Cov_.block(SG_ID, SG_ID, 3, 3) = imuerror_std.gyrscale.cwiseProduct(imuerror_std.gyrscale).asDiagonal();

    // 18，18，21，21：加速度计比例因子
    m_Cov_.block(SA_ID, SA_ID, 3, 3) = imuerror_std.accscale.cwiseProduct(imuerror_std.accscale).asDiagonal();
}

void GIEngine::newImuDataProcess()
{

    // 当前IMU时间作为系统当前状态时间,
    // set current IMU time as the current state time
    m_timeStamp_ = m_imuCur_.time;

    // 如果GNSS有效，则将更新时间设置为GNSS时间
    // set update time as the gnss time if gnssdata is valid
    double updateTime = m_gnssData_.isvalid ? m_gnssData_.time : -1;

    // 判断是否需要进行GNSS更新
    // determine if we should do GNSS update
    int res = isToUpdate(m_imuPre_.time, m_imuCur_.time, updateTime);

    if (res == 0)
    {
        // 更新时间不在imutimt1和imutime2之间，且不靠近任何一个
        // 只传播导航状态
        // only propagate navigation state
        insPropagation(m_imuPre_, m_imuCur_);
    }
    else if (res == 1)
    {
        // 更新时间靠近imutime1
        // GNSS数据靠近上一历元，先对上一历元进行GNSS更新
        // gnssdata is near to the previous imudata, we should firstly do gnss update
        gnssUpdate(m_gnssData_);
        stateFeedback();

        m_pvaPre_ = m_pvaCur_;
        insPropagation(m_imuPre_, m_imuCur_);
    }
    else if (res == 2)
    {
        // 更新时间靠近imutime2
        // GNSS数据靠近当前历元，先对当前IMU进行状态传播
        // gnssdata is near current imudata, we should firstly propagate navigation state
        insPropagation(m_imuPre_, m_imuCur_);
        gnssUpdate(m_gnssData_);
        stateFeedback();
    }
    else // res == 3
    {
        // 更新时间在imutime1和imutime2之间, 但不靠近任何一个
        // GNSS数据在两个IMU数据之间(不靠近任何一个), 将当前IMU内插到整秒时刻
        // gnssdata is between the two imudata, we interpolate current imudata to gnss time
        IMU midimu;
        imuInterpolate(m_imuPre_, m_imuCur_, updateTime, midimu);

        // 对前一半IMU进行状态传播
        // propagate navigation state for the first half imudata
        insPropagation(m_imuPre_, midimu);

        // 整秒时刻进行GNSS更新，并反馈系统状态
        // do GNSS position update at the whole second and feedback system states
        gnssUpdate(m_gnssData_);
        stateFeedback();

        // 对后一半IMU进行状态传播
        // propagate navigation state for the second half imudata
        m_pvaPre_ = m_pvaCur_;
        insPropagation(midimu, m_imuCur_);
    }

    // 检查协方差矩阵对角线元素
    // check diagonal elements of current covariance matrix
    checkCov();

    // 更新上一时刻的状态和IMU数据
    // update system state and imudata at the previous epoch
    m_pvaPre_ = m_pvaCur_;
    m_imuPre_ = m_imuCur_;
}

void GIEngine::imuCompensate(IMU& imu)
{

    // 补偿IMU零偏
    // compensate the imu bias
    imu.dtheta -= m_imuError_.gyrbias * imu.dt;
    imu.dvel -= m_imuError_.accbias * imu.dt;

    // 补偿IMU比例因子
    // compensate the imu scale
    Eigen::Vector3d gyrscale, accscale;
    gyrscale = Eigen::Vector3d::Ones() + m_imuError_.gyrscale;
    accscale = Eigen::Vector3d::Ones() + m_imuError_.accscale;
    imu.dtheta = imu.dtheta.cwiseProduct(gyrscale.cwiseInverse());
    imu.dvel = imu.dvel.cwiseProduct(accscale.cwiseInverse());
}

void GIEngine::insPropagation(IMU& imupre, IMU& imucur)
{

    // 对当前IMU数据(imucur)补偿误差, 上一IMU数据(imupre)已经补偿过了
    // compensate imu error to 'imucur', 'imupre' has been compensated
    imuCompensate(imucur);
    // IMU状态更新(机械编排算法)
    // update imustate(mechanization)
    INSMech::insMech(m_pvaPre_, m_pvaCur_, imupre, imucur);

    // 系统噪声传播，姿态误差采用phi角误差模型
    // system noise propagate, phi-angle error model for attitude error
    Eigen::MatrixXd Phi, F, Qd, G;

    // 初始化Phi阵(状态转移矩阵)，F阵，Qd阵(传播噪声阵)，G阵(噪声驱动阵)
    // initialize Phi (state transition), F matrix, Qd(propagation noise) and G(noise driven) matrix
    Phi.resizeLike(m_Cov_);
    F.resizeLike(m_Cov_);
    Qd.resizeLike(m_Cov_);
    G.resize(m_constRANK, m_constNOISERANK);
    Phi.setIdentity();
    F.setZero();
    Qd.setZero();
    G.setZero();

    // 使用上一历元状态计算状态转移矩阵
    // compute state transition matrix using the previous state
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    double gravity;
    rmrn = Earth::meridianPrimeVerticalRadius(m_pvaPre_.pos[0]);
    gravity = Earth::gravity(m_pvaPre_.pos);
    wie_n << WGS84_WIE * cos(m_pvaPre_.pos[0]), 0, -WGS84_WIE * sin(m_pvaPre_.pos[0]);
    wen_n << m_pvaPre_.vel[1] / (rmrn[1] + m_pvaPre_.pos[2]), -m_pvaPre_.vel[0] / (rmrn[0] + m_pvaPre_.pos[2]),
        -m_pvaPre_.vel[1] * tan(m_pvaPre_.pos[0]) / (rmrn[1] + m_pvaPre_.pos[2]);

    Eigen::Matrix3d temp;
    Eigen::Vector3d accel, omega;
    double rmh, rnh;

    rmh = rmrn[0] + m_pvaPre_.pos[2];
    rnh = rmrn[1] + m_pvaPre_.pos[2];
    accel = imucur.dvel / imucur.dt;
    omega = imucur.dtheta / imucur.dt;

    // 位置误差
    // position error
    temp.setZero();
    temp(0, 0) = -m_pvaPre_.vel[2] / rmh;
    temp(0, 2) = m_pvaPre_.vel[0] / rmh;
    temp(1, 0) = m_pvaPre_.vel[1] * tan(m_pvaPre_.pos[0]) / rnh;
    temp(1, 1) = -(m_pvaPre_.vel[2] + m_pvaPre_.vel[0] * tan(m_pvaPre_.pos[0])) / rnh;
    temp(1, 2) = m_pvaPre_.vel[1] / rnh;
    F.block(P_ID, P_ID, 3, 3) = temp;
    F.block(P_ID, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 速度误差
    // velocity error
    temp.setZero();
    temp(0, 0) = -2 * m_pvaPre_.vel[1] * WGS84_WIE * cos(m_pvaPre_.pos[0]) / rmh -
        pow(m_pvaPre_.vel[1], 2) / rmh / rnh / pow(cos(m_pvaPre_.pos[0]), 2);
    temp(0, 2) = m_pvaPre_.vel[0] * m_pvaPre_.vel[2] / rmh / rmh - pow(m_pvaPre_.vel[1], 2) * tan(m_pvaPre_.pos[0]) / rnh / rnh;
    temp(1, 0) = 2 * WGS84_WIE * (m_pvaPre_.vel[0] * cos(m_pvaPre_.pos[0]) - m_pvaPre_.vel[2] * sin(m_pvaPre_.pos[0])) / rmh +
        m_pvaPre_.vel[0] * m_pvaPre_.vel[1] / rmh / rnh / pow(cos(m_pvaPre_.pos[0]), 2);
    temp(1, 2) = (m_pvaPre_.vel[1] * m_pvaPre_.vel[2] + m_pvaPre_.vel[0] * m_pvaPre_.vel[1] * tan(m_pvaPre_.pos[0])) / rnh / rnh;
    temp(2, 0) = 2 * WGS84_WIE * m_pvaPre_.vel[1] * sin(m_pvaPre_.pos[0]) / rmh;
    temp(2, 2) = -pow(m_pvaPre_.vel[1], 2) / rnh / rnh - pow(m_pvaPre_.vel[0], 2) / rmh / rmh +
        2 * gravity / (sqrt(rmrn[0] * rmrn[1]) + m_pvaPre_.pos[2]);
    F.block(V_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 0) = m_pvaPre_.vel[2] / rmh;
    temp(0, 1) = -2 * (WGS84_WIE * sin(m_pvaPre_.pos[0]) + m_pvaPre_.vel[1] * tan(m_pvaPre_.pos[0]) / rnh);
    temp(0, 2) = m_pvaPre_.vel[0] / rmh;
    temp(1, 0) = 2 * WGS84_WIE * sin(m_pvaPre_.pos[0]) + m_pvaPre_.vel[1] * tan(m_pvaPre_.pos[0]) / rnh;
    temp(1, 1) = (m_pvaPre_.vel[2] + m_pvaPre_.vel[0] * tan(m_pvaPre_.pos[0])) / rnh;
    temp(1, 2) = 2 * WGS84_WIE * cos(m_pvaPre_.pos[0]) + m_pvaPre_.vel[1] / rnh;
    temp(2, 0) = -2 * m_pvaPre_.vel[0] / rmh;
    temp(2, 1) = -2 * (WGS84_WIE * cos(m_pvaPre_.pos(0)) + m_pvaPre_.vel[1] / rnh);
    F.block(V_ID, V_ID, 3, 3) = temp;
    F.block(V_ID, PHI_ID, 3, 3) = Rotation::skewSymmetric(m_pvaPre_.att.cbn * accel);
    F.block(V_ID, BA_ID, 3, 3) = m_pvaPre_.att.cbn;
    F.block(V_ID, SA_ID, 3, 3) = m_pvaPre_.att.cbn * (accel.asDiagonal());

    // 姿态误差
    // attitude error
    temp.setZero();
    temp(0, 0) = -WGS84_WIE * sin(m_pvaPre_.pos[0]) / rmh;
    temp(0, 2) = m_pvaPre_.vel[1] / rnh / rnh;
    temp(1, 2) = -m_pvaPre_.vel[0] / rmh / rmh;
    temp(2, 0) = -WGS84_WIE * cos(m_pvaPre_.pos[0]) / rmh - m_pvaPre_.vel[1] / rmh / rnh / pow(cos(m_pvaPre_.pos[0]), 2);
    temp(2, 2) = -m_pvaPre_.vel[1] * tan(m_pvaPre_.pos[0]) / rnh / rnh;
    F.block(PHI_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 1) = 1 / rnh;
    temp(1, 0) = -1 / rmh;
    temp(2, 1) = -tan(m_pvaPre_.pos[0]) / rnh;
    F.block(PHI_ID, V_ID, 3, 3) = temp;
    F.block(PHI_ID, PHI_ID, 3, 3) = -Rotation::skewSymmetric(wie_n + wen_n);
    F.block(PHI_ID, BG_ID, 3, 3) = -m_pvaPre_.att.cbn;
    F.block(PHI_ID, SG_ID, 3, 3) = -m_pvaPre_.att.cbn * (omega.asDiagonal());

    // IMU零偏误差和比例因子误差，建模成一阶高斯-马尔科夫过程
    // imu bias error and scale error, modeled as the first-order Gauss-Markov process
    F.block(BG_ID, BG_ID, 3, 3) = -1 / m_options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(BA_ID, BA_ID, 3, 3) = -1 / m_options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SG_ID, SG_ID, 3, 3) = -1 / m_options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SA_ID, SA_ID, 3, 3) = -1 / m_options_.imunoise.corr_time * Eigen::Matrix3d::Identity();

    // 系统噪声驱动矩阵
    // system noise driven matrix
    G.block(V_ID, VRW_ID, 3, 3) = m_pvaPre_.att.cbn;
    G.block(PHI_ID, ARW_ID, 3, 3) = m_pvaPre_.att.cbn;
    G.block(BG_ID, BGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(BA_ID, BASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SG_ID, SGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SA_ID, SASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 状态转移矩阵
    // compute the state transition matrix
    Phi.setIdentity();
    Phi = Phi + F * imucur.dt;

    // 计算系统传播噪声
    // compute system propagation noise
    Qd = G * m_Qc_ * G.transpose() * imucur.dt;
    Qd = (Phi * Qd * Phi.transpose() + Qd) / 2;

    // EKF预测传播系统协方差和系统误差状态
    // do EKF predict to propagate covariance and error state
    EKFPredict(Phi, Qd);
}

void GIEngine::gnssUpdate(GNSS& gnssdata)
{

    // IMU位置转到GNSS天线相位中心位置
    // convert IMU position to GNSS antenna phase center position
    Eigen::Vector3d antenna_pos;
    Eigen::Matrix3d Dr, Dr_inv;
    Dr_inv = Earth::DRi(m_pvaCur_.pos);
    Dr = Earth::DR(m_pvaCur_.pos);
    antenna_pos = m_pvaCur_.pos + Dr_inv * m_pvaCur_.att.cbn * m_options_.antlever;

    // GNSS位置测量新息
    // compute GNSS position innovation
    Eigen::MatrixXd dz;
    dz = Dr * (antenna_pos - gnssdata.blh);

    // 构造GNSS位置观测矩阵
    // construct GNSS position measurement matrix
    Eigen::MatrixXd H_gnsspos;
    H_gnsspos.resize(3, m_Cov_.rows());
    H_gnsspos.setZero();
    H_gnsspos.block(0, P_ID, 3, 3) = Eigen::Matrix3d::Identity();
    H_gnsspos.block(0, PHI_ID, 3, 3) = Rotation::skewSymmetric(m_pvaCur_.att.cbn * m_options_.antlever);

    // 位置观测噪声阵
    // construct measurement noise matrix
    Eigen::MatrixXd R_gnsspos;
    R_gnsspos = gnssdata.std.cwiseProduct(gnssdata.std).asDiagonal();

    // EKF更新协方差和误差状态
    // do EKF update to update covariance and error state
    EKFUpdate(dz, H_gnsspos, R_gnsspos);

    // GNSS更新之后设置为不可用
    // Set GNSS invalid after update
    gnssdata.isvalid = false;
}

int GIEngine::isToUpdate(double imutime1, double imutime2, double updatetime) const
{

    if (abs(imutime1 - updatetime) < m_constTIME_ALIGN_ERR)
    {
        // 更新时间靠近imutime1
        // updatetime is near to imutime1
        return 1;
    }
    else if (abs(imutime2 - updatetime) <= m_constTIME_ALIGN_ERR)
    {
        // 更新时间靠近imutime2
        // updatetime is near to imutime2
        return 2;
    }
    else if (imutime1 < updatetime && updatetime < imutime2)
    {
        // 更新时间在imutime1和imutime2之间, 但不靠近任何一个
        // updatetime is between imutime1 and imutime2, but not near to either
        return 3;
    }
    else
    {
        // 更新时间不在imutimt1和imutime2之间，且不靠近任何一个
        // updatetime is not bewteen imutime1 and imutime2, and not near to either.
        return 0;
    }
}

void GIEngine::EKFPredict(Eigen::MatrixXd& Phi, Eigen::MatrixXd& Qd)
{

    assert(Phi.rows() == m_Cov_.rows());
    assert(Qd.rows() == m_Cov_.rows());

    // 传播系统协方差和误差状态
    // propagate system covariance and error state
    m_Cov_ = Phi * m_Cov_ * Phi.transpose() + Qd;
    m_dx_ = Phi * m_dx_;
}

void GIEngine::EKFUpdate(Eigen::MatrixXd& dz, Eigen::MatrixXd& H, Eigen::MatrixXd& R)
{

    assert(H.cols() == m_Cov_.rows());
    assert(dz.rows() == H.rows());
    assert(dz.rows() == R.rows());
    assert(dz.cols() == 1);

    // 计算Kalman增益
    // Compute Kalman Gain
    auto temp = H * m_Cov_ * H.transpose() + R;
    Eigen::MatrixXd K = m_Cov_ * H.transpose() * temp.inverse();

    // 更新系统误差状态和协方差
    // update system error state and covariance
    Eigen::MatrixXd I;
    I.resizeLike(m_Cov_);
    I.setIdentity();
    I = I - K * H;
    // 如果每次更新后都进行状态反馈，则更新前dx_一直为0，下式可以简化为：dx_ = K * dz;
    // if state feedback is performed after every update, dx_ is always zero before the update
    // the following formula can be simplified as : dx_ = K * dz;
    m_dx_ = m_dx_ + K * (dz - H * m_dx_);
    m_Cov_ = I * m_Cov_ * I.transpose() + K * R * K.transpose();
}

void GIEngine::stateFeedback()
{

    Eigen::Vector3d vectemp;

    // 位置误差反馈
    // posisiton error feedback
    Eigen::Vector3d delta_r = m_dx_.block(P_ID, 0, 3, 1);
    Eigen::Matrix3d Dr_inv = Earth::DRi(m_pvaCur_.pos);
    m_pvaCur_.pos -= Dr_inv * delta_r;

    // 速度误差反馈
    // velocity error feedback
    vectemp = m_dx_.block(V_ID, 0, 3, 1);
    m_pvaCur_.vel -= vectemp;

    // 姿态误差反馈
    // attitude error feedback
    vectemp = m_dx_.block(PHI_ID, 0, 3, 1);
    Eigen::Quaterniond qpn = Rotation::rotvec2quaternion(vectemp);
    m_pvaCur_.att.qbn = qpn * m_pvaCur_.att.qbn;
    m_pvaCur_.att.cbn = Rotation::quaternion2matrix(m_pvaCur_.att.qbn);
    m_pvaCur_.att.euler = Rotation::matrix2euler(m_pvaCur_.att.cbn);

    // IMU零偏误差反馈
    // IMU bias error feedback
    vectemp = m_dx_.block(BG_ID, 0, 3, 1);
    m_imuError_.gyrbias += vectemp;
    vectemp = m_dx_.block(BA_ID, 0, 3, 1);
    m_imuError_.accbias += vectemp;

    // IMU比例因子误差反馈
    // IMU sacle error feedback
    vectemp = m_dx_.block(SG_ID, 0, 3, 1);
    m_imuError_.gyrscale += vectemp;
    vectemp = m_dx_.block(SA_ID, 0, 3, 1);
    m_imuError_.accscale += vectemp;

    // 误差状态反馈到系统状态后,将误差状态清零
    // set 'dx' to zero after feedback error state to system state
    m_dx_.setZero();
}

NavState GIEngine::getNavState()
{

    NavState state;

    state.pos = m_pvaCur_.pos;
    state.vel = m_pvaCur_.vel;
    state.euler = m_pvaCur_.att.euler;
    state.imuerror = m_imuError_;

    return state;
}
