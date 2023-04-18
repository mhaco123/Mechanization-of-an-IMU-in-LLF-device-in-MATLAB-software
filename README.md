# Mechanization-of-an-IMU-in-LLF-device-in-MATLAB-software
Mechanization is the process of converting IMU output into position, speed and direction information. The outputs include the rotational speed about the three body axes b measured by the triple gyroscope and the three specific forces f along the body axes measured by the triple accelerometer.
  Mechanization is a recursive process that starts with a specific set of initial values and iterates over the output.
  
  The first and second step: measure the primary data and reduce the noise from the raw data
  The third step: calculating and updating the period matrix
  Fourth step: Orientation calculation
  Step five: Calculate the speed
  The sixth step: calculate the position
  
  Evaluation and validation
To check the results, four scenarios designed for verification have been used. In this part, using the Logger Sensor Android software, we extract the real data of the accelerometer and GPS in the form of a csv format file and put them in the relations.

First scenario: a stationary object at a certain height
We expect it to be able to show changes in the movement of the earth; Due to the change of latitude and longitude, the output equations must change in line with the rotation of the earth. That is, we assume that the height changes are equal to zero

The second scenario: the fall of a stationary object at a certain height
In the next test from a height of 1000 meters and based on the data taken from the Android software, we expect it to reach the ground after a time of 4.14.
Assuming a height of 160 meters, we expect the height to reach zero within a certain time based on the law of gravity

The third scenario: moving a straight line on the ground based on real data. In this scenario, we expect the longitude and latitude to change according to the given data.
