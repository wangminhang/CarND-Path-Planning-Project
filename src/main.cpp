#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// Param to be adjusted
const double TARGET_SPEED = 49.5; // max speed
const double TARGET_ACC = 6;	// max acceleration
const double TARGET_JERK = 4.5;   // max jerk
const double REFRESH_TIME = 0.02;		// refresh time between two position
const double LANE_WIDTH = 4;			// width of lane

// For keeping track externally of lane change
bool g_changing_lane = false;
int g_target_lane = 0;
double time_horizon = 1;
double critical_distance = 6;					// distance to slow down
double prepare_lane_change = 30;				// anticipate by trying to change lane early
double safety_lane_change_distance_back = 8;   // distance from back car when want change line
double safety_lane_change_distance_front = 15; // distance from front car when want change line

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y) {
    int prev_wp = -1;

    while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
        prev_wp++;
    }

    int wp2 = (prev_wp + 1) % maps_x.size();

    double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s - maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
    double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

    double perp_heading = heading - pi() / 2;

    double x = seg_x + d * cos(perp_heading);
    double y = seg_y + d * sin(perp_heading);

    return {x, y};
}

// helper functions about coordinate transform
vector<double> global_to_local(const double &x, const double &y, const double &x_ref, const double &y_ref, const double &theta_ref)
{
	double xl, yl, xc, yc;

	xc = x - x_ref;
	yc = y - y_ref;

	xl = xc * cos(deg2rad(theta_ref)) + yc * sin(deg2rad(theta_ref));
	yl = yc * cos(deg2rad(theta_ref)) - xc * sin(deg2rad(theta_ref));

	return {xl, yl};
}

void global_to_local(vector<double> &x, vector<double> &y, const double &x_ref, const double &y_ref, const double &theta_ref)
{
	for (int i = 0; i < x.size(); ++i)
	{
		auto z = global_to_local(x[i], y[i], x_ref, y_ref, theta_ref);
		x[i] = z[0];
		y[i] = z[1];
	}
}

vector<double> local_to_global(const double &x, const double &y, const double &x_ref, const double &y_ref, const double &theta_ref)
{
	double xg, yg;

	xg = x * cos(deg2rad(theta_ref)) - y * sin(deg2rad(theta_ref));
	yg = y * cos(deg2rad(theta_ref)) + x * sin(deg2rad(theta_ref));

	xg += x_ref;
	yg += y_ref;

	return {xg, yg};
}

void local_to_global(vector<double> &x, vector<double> &y, const double &x_ref, const double &y_ref, const double &theta_ref)
{
	for (int i = 0; i < x.size(); ++i)
	{
		auto z = local_to_global(x[i], y[i], x_ref, y_ref, theta_ref);
		x[i] = z[0];
		y[i] = z[1];
	}
}

void check_close_and_line_change(vector<vector<double>> sensor_fusion, double car_s, double my_future_car_s, bool &too_close, bool &try_lane_change)
{
	for (auto &other_car : sensor_fusion)
	{
		double other_car_lane = (other_car[6] + LANE_WIDTH / 2) / LANE_WIDTH;
		if (abs(other_car_lane - g_target_lane) < 0.7) // in the same lane or getting close
		{
			double other_car_s = other_car[5];
			double other_car_speed = sqrt(other_car[3] * other_car[3] + other_car[4] * other_car[4]);
			double other_future_car_s = other_car_s + other_car_speed * 1609 / 3600 * time_horizon;
			if ((((other_future_car_s - critical_distance) < my_future_car_s) || ((other_car_s - critical_distance) < car_s)) && (car_s < other_car_s)) // we are either already too close or about to be too close
			{
				too_close = true;
				try_lane_change = true;
				break;
			}
			if ((((other_future_car_s - prepare_lane_change) < my_future_car_s) || ((other_car_s - prepare_lane_change) < car_s)) && (car_s < other_car_s)) // we are either already too close or about to be too close
			{
				try_lane_change = true;
			}
		}
	}
}

void check_other_line_safe(vector<vector<double>> sensor_fusion, double car_s, double my_future_car_s, bool &left_line_safe, bool &right_line_safe)
{
	for (auto &one_car : sensor_fusion)
	{
		double one_car_s = one_car[5];
		double one_car_speed = sqrt(one_car[3] * one_car[3] + one_car[4] * one_car[4]);
		double other_future_car_s = one_car_s + one_car_speed * 1609 / 3600 * time_horizon;

		// if current one car will override my car, then can not change line
		if ((((one_car_s - car_s) > -safety_lane_change_distance_back) && ((one_car_s - car_s) < safety_lane_change_distance_front)) ||
			(((other_future_car_s - my_future_car_s) > -safety_lane_change_distance_back) && ((other_future_car_s - my_future_car_s) < safety_lane_change_distance_front)) ||
			((one_car_s < car_s - safety_lane_change_distance_back) && (other_future_car_s > my_future_car_s + safety_lane_change_distance_front)) ||
			((other_future_car_s < my_future_car_s - safety_lane_change_distance_back) && (one_car_s > car_s + safety_lane_change_distance_front)))
		{
			double one_car_lane = (one_car[6] + LANE_WIDTH / 2) / LANE_WIDTH;
			if (((one_car_lane - g_target_lane) < -0.5) && ((one_car_lane - g_target_lane) > -1.5)) // other car is on our left
			{
				left_line_safe = false;
			}
			if (((one_car_lane - g_target_lane) > 0.5) && ((one_car_lane - g_target_lane) < 1.5)) // other car is on our right
			{
				right_line_safe = false;
			}
		}
	}
}

// generate line
vector<vector<double>> line_gen(const double &car_s, const double &car_d, const double &car_x, const double &car_y, const double &car_yaw, const double &car_speed,
								   const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y, const vector<double> &maps_dx, const vector<double> &maps_dy,
								   const vector<double> &previous_path_x, const vector<double> &previous_path_y,
								   vector<vector<double>> sensor_fusion)
{
	vector<double> pos_x, pos_y;

	int keep_points = 10;			// Number of points from previous trajectory we will keep
	int spline_prev_trajectory = 3; // Number of points from previous trajectory we use for spline
	int previous_path_size = previous_path_x.size();
	int prev_points_used = min(previous_path_size, keep_points);

	// Define target lane if it has not been initialized yet and save it in global variable
	// will use current line as target lane
	if (!g_target_lane)
	{
		g_target_lane = int(round((car_d + LANE_WIDTH / 2.0) / LANE_WIDTH));
	}
	double target_d = LANE_WIDTH * (g_target_lane - 1 / 2.0);
	double target_horizon = 15.;

	// See if we will be too close to other cars in 3 seconds based on current velocities
	bool too_close = false;
	bool try_change_lane = false;
	double my_future_car_s = car_s + car_speed * 1609 / 3600 * time_horizon;
	check_close_and_line_change(sensor_fusion, car_s, my_future_car_s, too_close, try_change_lane);

	// Check if we want to change lane as long as we are not already changing lane
	if (try_change_lane && !g_changing_lane)
	{
		bool left_lane_free = true;
		bool right_lane_free = true;
		if (g_target_lane == 1)
		{
			left_lane_free = false;
		}
		if (g_target_lane == 3)
		{
			right_lane_free = false;
		}

		// scan all the cars from sensor fusion
		check_other_line_safe(sensor_fusion, car_s, my_future_car_s, left_lane_free, right_lane_free);

		// Change lane if one is free
		if (left_lane_free)
		{
			g_changing_lane = true;
			g_target_lane -= 1;
		}
		else if (right_lane_free)
		{
			g_changing_lane = true;
			g_target_lane += 1;
		}
	}

	if (g_changing_lane)
	{
		target_horizon = 35.0; // This ensure we don't change lane too fast
		target_d = LANE_WIDTH * (g_target_lane - 1.0 / 2.0);

		// Check if we are done changing lane
		if (abs((car_d + LANE_WIDTH / 2.0) / LANE_WIDTH - g_target_lane) < 0.25)
		{
			g_changing_lane = false; // we have finished our change of lane
		}
	}

	// Define our target speed
	double target_speed = car_speed * 1609 / 3600;
	if (too_close)
	{
		target_speed -= TARGET_ACC * REFRESH_TIME * keep_points;
	}
	else
	{
		target_speed += TARGET_ACC * REFRESH_TIME * keep_points;
		target_speed = min(target_speed, TARGET_SPEED * 1609 / 3600);
	}

	// Create a spline based on 3 points: current position and 2 points front and back of the car
	tk::spline s_xy;

	if (previous_path_size > spline_prev_trajectory)
	{
		std::vector<double> s_X(spline_prev_trajectory + 2), s_Y(spline_prev_trajectory + 2);

		// add last trajectory points
		for (int i = 0; i < spline_prev_trajectory; ++i)
		{
			s_X[i] = previous_path_x[prev_points_used - spline_prev_trajectory + i];
			s_Y[i] = previous_path_y[prev_points_used - spline_prev_trajectory + i];
		}

		// add points based on last point from target trajectory + a small distance (target_horizon)
		double last_x = s_X[spline_prev_trajectory - 1];
		double last_y = s_Y[spline_prev_trajectory - 1];
		double last_x2 = s_X[spline_prev_trajectory - 2];
		double last_y2 = s_Y[spline_prev_trajectory - 2];
		double last_theta = atan2(last_y - last_y2, last_x - last_x2);
		auto frenet = getFrenet(last_x, last_y, last_theta, maps_x, maps_y);
		auto next_z = getXY(frenet[0] + target_horizon, target_d, maps_s, maps_x, maps_y);
		s_X[spline_prev_trajectory] = next_z[0];
		s_Y[spline_prev_trajectory] = next_z[1];
		next_z = getXY(frenet[0] + 2 * target_horizon, target_d, maps_s, maps_x, maps_y);
		s_X[spline_prev_trajectory + 1] = next_z[0];
		s_Y[spline_prev_trajectory + 1] = next_z[1];

		// Convert all points to local coordinates
		global_to_local(s_X, s_Y, car_x, car_y, car_yaw);

		// create the spline
		s_xy.set_points(s_X, s_Y);
	}
	else
	// We don't have enough points so we will use last waypoint, current position and next waypoint
	{
		std::vector<double> s_X(3), s_Y(3);
		auto prev_z = getXY(car_s - target_horizon, target_d, maps_s, maps_x, maps_y);
		auto next_z = getXY(car_s + target_horizon, target_d, maps_s, maps_x, maps_y);
		s_X[0] = prev_z[0];
		s_Y[0] = prev_z[1];
		s_X[1] = car_x;
		s_Y[1] = car_y;
		s_X[2] = next_z[0];
		s_Y[2] = next_z[1];

		// Convert all points to local coordinates
		global_to_local(s_X, s_Y, car_x, car_y, car_yaw);

		// create the spline
		s_xy.set_points(s_X, s_Y);
	}

	// Keep part of the previous trajectory if possible
	for (int i = 0; i < prev_points_used; ++i)
	{
		pos_x.push_back(previous_path_x[i]);
		pos_y.push_back(previous_path_y[i]);
	}

	// Get last point we took from previous trajectory, use current car position if not available
	double last_x = car_x;
	double last_y = car_y;
	if (prev_points_used)
	{
		last_x = pos_x.back();
		last_y = pos_y.back();
	}
	auto last_z_local = global_to_local(last_x, last_y, car_x, car_y, car_yaw);
	double last_x_local = last_z_local[0];
	double last_y_local = last_z_local[1];

	// Add new points based on target speed
	double distance_inc = target_speed * REFRESH_TIME; // target distance in m
	for (int i = pos_x.size(); i < 50; ++i)			   // have a total of 50 elements, equivalent to 1 sec
	{
		// Look for (x, y) approximately at distance_inc
		double next_x_local = last_x_local + distance_inc;
		double next_y_local = s_xy(next_x_local);
		double d = distance(last_x_local, last_y_local, next_x_local, next_y_local);

		// Readjust (x, y) based on calculated distance
		next_x_local = last_x_local + (next_x_local - last_x_local) * distance_inc / d;
		next_y_local = s_xy(next_x_local);

		// Add calculated positions to trajectory
		auto next_z_global = local_to_global(next_x_local, next_y_local, car_x, car_y, car_yaw);
		pos_x.push_back(next_z_global[0]);
		pos_y.push_back(next_z_global[1]);

		// Update last points of current trajectory
		last_x_local = next_x_local;
		last_y_local = next_y_local;
	}

	return {pos_x, pos_y};
}

int main()
{
	uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

			auto generated_path = line_gen(car_s, car_d, car_x, car_y, car_yaw, car_speed, map_waypoints_s, map_waypoints_x, map_waypoints_y, map_waypoints_dx, map_waypoints_dy, previous_path_x, previous_path_y, sensor_fusion);
			next_x_vals = generated_path[0];
			next_y_vals = generated_path[1];

          	// define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
