from time import time
from math import sqrt, pi, sin, cos, atan2
from numpy import array, dot, zeros
from numpy.linalg import norm

from Robot3WD.PID import PIDControllers

from myRobot import *

# Set up sampling period T_s and stopping time T_f
T_s = 0.02 
T_f = 20.0

# Set up initial pose
x_0 = 0.0
y_0 = 0.0
theta_0 = (np.pi) / 6

# Set up goal position
x_f = 0.7
y_f = 1.0
theta_f = (np.pi)/2

# Set up p_0 and p_f
p_0 = array([x_0, y_0, theta_0]).T
p_f = array([x_f, y_f, theta_f]).T

# Set up error tolerance
epsilon = sqrt( 2.0*(0.5/100.0)**2)

# set up d_free and d_min
d_min = 0.08
d_free = 1.25*d_min

# Set up controller gains
k_rho = 0.5
k_beta = -0.5
k_alpha = (2.0/pi)*k_rho - (5.0/3.0)*k_beta + 0.5

#Set up PID

Kp = array([0.2, 0.2, 0.2]).T
Ki = array([0.2, 0.2, 0.2]).T
Kd = array([0.01, 0.01, 0.01]).T
pid = PIDControllers(Kp, Ki, Kd, T_s)

# Initialize vector pr_dot to be used for inverse kinematics
pr_dot = array([0, 0, 0]).T

# Initialize vector d for storing the sensor measurements d_i 
d = array([10.0, 10.0, 10.0, 10.0, 10.0, 10.0]).T

# Initialize a 6x3 matrix u_iR.  The vector u_i^R is assigned the ith column
# of this matrix.  See Eq. (4.6) in lab manual.
u_iR = zeros(shape=(3,6))

# Sensor frame location
R = 0.13
halfR = R/2.0
root3Rover2 = sqrt(3.0)*R/2.0

# Sensor frame location
sensor_loc = np.array([
   [-halfR, -root3Rover2, -(2.0/3.0)*pi],
   [-root3Rover2,halfR, (5.0/6.0)*pi],
   [ 0,R, 0.5*pi],
   [ halfR,root3Rover2,(1.0/3.0)*pi],
   [root3Rover2,halfR,(1.0/6.0)*pi ],
   [root3Rover2,-halfR,(-1.0/6.0)*pi]
])

# open a file to store robot pose for plotting later
f = open('pose.csv', 'w')

# Initial robot and other variables
robot = myRobot(T_s)

robot.initialize(theta_0)

p = p_0
dp = p_f - p_0
rho = norm(dp[0:2])
goal_reached = False
elapsed_time = 0.0

# Set start time
start_time = time()
neverdone = True
# Control loop
while ( (not goal_reached) and (elapsed_time < T_f) ):

    robot.get_readings_update()

    # save data to file for plotting
    f.write('%10.3e %10.3e % 10.3e %10.3e %10.3e\n' % (elapsed_time, p[0], p[1], p[2], rho))

    theta = robot.orientation
    
    # Use forward kinematics to get p_dot
    p_dot = robot.forward_kinematics(robot.angular_velocities, theta)
    pd_dot = pid(p_f, p)
    
    # Get sensor measurements w.r.t. sensor frames
    for i in range(0,6):
       d[i] = robot.Vraw_to_distance(robot.ir_sensors_raw_values[i])

    # Check if obstacle is present
    # If true, determine temporary goal position and use it to calculate
    # rho bar and alpha bar.  See Eqs. 4.6 to 4.11.
    d_T = 0.0
    HR0 = robot.HMatrix(p)
    p = (p + p_dot*T_s)
    u_r = array([0.0,0.0,0.0]).T
    rho = sqrt(((p_f[0]-p[0])**2)+((p_f[1]-p[1])**2)) 
    alpha = atan2((p_f[1]-p[1]), (p_f[0]-p[0])) - robot.orientation

    if min(d) < d_min:
        d_T = sum(i for i in d if i > d_free)
        gen = (x for x in d if x > d_free)
        for x in gen:
            sidx = np.where(d == x)
            sidx = np.squeeze(sidx)
            if sidx.size > 1 :
                d[sidx[0]] = 0.0
                sidx = np.delete(sidx, 1)
                sid = sidx[0]
            else:
                sid = sidx
            print("S",sid,": ",x)
            u_i = sensor_loc[sid,:]/ np.linalg.norm(sensor_loc[sid,:])
            u_i = np.squeeze(u_i)
            fac = x/d_T
            u_r+= np.dot(HR0, fac*u_i)
        u_0 = u_r
        p_tmp = d_free*u_0
        
        p_bar = sqrt(((p_tmp[0]-p[0])**2)+((p_tmp[1]-p[1])**2))
        a_bar = atan2((p_tmp[1]-p[1]), (p_tmp[0]-p[0])) -theta 
        alpha = a_bar
        rho = p_bar

    beta = -(theta + alpha)

    # Determine linear and angular velocities of the robot body (v, w)
    # using the control law given in Eqs. (4.12) and (4.13)
    v = k_rho*rho
    w = k_alpha*alpha + k_beta*beta

    # Determine pr_dot
    pr_dot = array([v*cos(theta), v*sin(theta), w]).T
    # Now use Inverse Kinematics to determine wheel ref velocities
    wheel_ref_vel = robot.inverse_kinematics(pr_dot, theta)
    # Apply motor limits
    wheel_ref_vel = robot.motor_sat(wheel_ref_vel, (5 * np.pi))
    # Execute motion
    robot.set_angular_velocities(wheel_ref_vel)

    # Odometry update
    p = (p + p_dot*T_s)
    print("pos: ",p)

    # Replace calculated update for theta with measured value from IMU
    p[2]= robot.orientation

    # Check to see if goal is reached
    dp = p_f - p
    trho = norm(dp[0:2])
    goal_reached = ( trho <= epsilon)

    # time update
    elapsed_time = time() - start_time


print('ELPASED TIME = %s rho = %s' % (elapsed_time, rho))

# Either goal is reached or current_time > T_f
if goal_reached:
   print('Goal is reached, error norm is', rho)
else:
   print('Failed to reach goal, error norm is', rho)
print('Final position:', p)
print('Desired position:',p_f)
robot.stop()
robot.close()

f.close()
