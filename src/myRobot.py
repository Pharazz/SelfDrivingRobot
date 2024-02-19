import numpy as np

from math import sin, cos, exp, sqrt, pi
from numpy import array, dot

from Robot3WD.Robot3WD import Robot

class myRobot(Robot):
    def __init__(self, sampling_period, wheel_radius=None, L=None):
        Robot.__init__(self, sampling_period, wheel_radius, L)

    def forward_kinematics(self, wheel_angular_velocities, theta):
       L = self._L
       wheel_radius = self._wheel_radius
       print(wheel_radius)
       trig = np.array([[(2*np.sin(theta)),(2*np.cos(theta + ((np.pi)/6))), np.negative(2*np.sin(theta + ((np.pi)/3)))],
                          [np.negative(2*np.cos(theta)), 2 * np.cos(theta - ((np.pi)/3)), (2*np.cos(theta + ((np.pi)/3)))],
                          [np.negative(1/L), np.negative(1/L), np.negative(1/L)]])
       trig = trig*(wheel_radius/3)
       p_dot = np.dot(trig, wheel_angular_velocities)
       print(p_dot)
       return p_dot
        
    def inverse_kinematics(self, p_dot, theta):
        L = self._L
        wheel_radius = self._wheel_radius
        trig = np.array([[np.sin(theta), np.negative((np.cos(theta))), np.negative(L)],
                          [np.cos((np.pi) / 6 + theta), np.sin((np.pi) / 6 + theta), np.negative(L)],
                          [np.negative(np.cos((np.pi) / 6 - theta)), np.sin((np.pi) / 6 - theta), np.negative(L)]])
        wheel_angular_velocities = np.dot(trig, p_dot)
        wheel_angular_velocities = wheel_angular_velocities/wheel_radius
        return wheel_angular_velocities

    def move_left(self, vx, theta):
        p_dot = array([-vx, 0.0, 0.0]).T
        w = self.inverse_kinematics(p_dot, theta)
        self.set_angular_velocities(w)

    def move_forward(self, vy, theta):
        
         p_dot = array([0.0, vy, 0.0]).T
         w = self.inverse_kinematics(p_dot, theta)
         #print(type(w))
         self.set_angular_velocities(w)
        
    def move_backward(self, vy, theta):

         p_dot = array([0.0, -vy, 0.0]).T
         w = self.inverse_kinematics(p_dot, theta)
         self.set_angular_velocities(w)
        
    def move_right(self, vx, theta):
        
         p_dot = array([vx, 0.0, 0.0]).T
         w = self.inverse_kinematics(p_dot, theta)
         self.set_angular_velocities(w)
        
        
    def rotate_CCW(self, w, theta):
       p_dot = array([0.0, 0.0, w]).T
       w = self.inverse_kinematics(p_dot, theta)
       self.set_angular_velocities(w)
        
    def rotate_CW(self, w, theta):
        p_dot = array([0.0, 0.0, -w]).T
        w = self.inverse_kinematics(p_dot, theta)
        self.set_angular_velocities(w)

    def motor_sat(self, wheel_angular_velocities, limit_value):	
       wheel_angular_velocities[(wheel_angular_velocities < (-1.0*limit_value))] = (-1.0*limit_value)     
       wheel_angular_velocities[(wheel_angular_velocities > (limit_value))] = (limit_value)     
       return wheel_angular_velocities

    def HMatrix(self, q):
        H = array([[ cos(q[2]), -sin(q[2]),q[0] ], [sin(q[2]) ,cos(q[2]) ,q[1] ], [0 ,0 ,1.0 ]])
        return H

    def Vraw_to_distance(self, Vraw):
        c2 = 0.053956773
        c1 =  0.620876659
        d = c1*exp(-c2*sqrt(Vraw))
        return d

#### end of myrobot Class ###

#### end of myrobot.py ###
