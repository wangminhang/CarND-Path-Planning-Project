# Udacity Self-Driving Car Engineer Nanodegree
## Path Planning Project

## Infomation

The goal of this project was to implement a path planning algorithm that utilizes given way points and nearby vehicles (simulating sensor fusion).  The algorithm should enable the ego vehicle to drive around the track safely, in that the ego vehicle:

- Shall not collide with any other vehicle
- Shall not go off the track (3 far right lanes)
- Shall be able to overtake when the adjacent lane is free and the current lane has slower moving traffic
- Shall not exert any extreme accelerations or jerks in both lateral and longitudinal directions 

### Algorithm

The algorithm includes several sub modules:

- Construct interpolated waypoints of nearby area  using spline interpolation
- Determine ego car parameters and construct vehicle object
- Generate predictions of each vehicles and determining the next possible states of the ego vehicle
- Calculating the cost of each potential next state and selecting the state with the lowest cost
- Creating and using the trajectory of the lowest state so that the ego vehicle proceeds in the allocated direction and speed

### State
State include keep line, change lane and currently changing line.We keep track of the state with two variables:
* *g_changing_lane* will be true if we have decided to change lane and go back to false when our change of lane is complete,
* *g_target_lane* will represent which lane the car will follow, which is useful to generate consistent trajectories.

### Trajectory
For car's trajectory, we need think about two thing: should it change and how it change. When we decide the action, we should create a *spline* to generate the car's trajetory. Through these action, we should follow these principles:
* At each step, we keep the last 10 time frames from previous calculated trajectory to ensure safety of our trajectory, which is equivalent to 0.2 sec.
* At the begin of decision, we alway calculate the other car to check those car is far away enough from our car.
* We use a *simple cubic spline interpolation library* to generate the car's trajetory
* The target trajectory and the target speed have a prediction over a total of 50 time frames (1 second).

### Speed
The speed should follow these principles:
* Can not over max speed limit(50mph).
* Through a checking function to decide speed up or slow down. We add up or slow the speed by limit acceleration in the function. The logic is below:
  * go through all related car from the sensor fusion.
  * calculate the related value of those cars.
  * if the line of the other car is in our line, we check if they are in front of us and close.
  * we check our car is still with in a safe distance at future.
  * if the case is yes, we speed up our car, or slow down.

### Lane changing
To change lane safely, we should follow these steps :
* Check the safety of the lines we can change to.
* Compare the other cars from fusion data with our car and calculate the line we are in, then calculate which line we can change to.
* The comparison used is deciding whether coincident may happen. If the coincidence is occurred, we can not change line.
* The state of change lane is continued action. We keep the state by a variable.