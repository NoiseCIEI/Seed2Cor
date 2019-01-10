#ifndef MAP_H
#define MAP_H

#include "Array2D.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>

#ifndef POINT
#define POINT
template < class T >
class Point {
protected:
   T lon, lat;
public:
   Point(T lonin = -12345., T latin = -12345.) 
      : lon(lonin), lat(latin) {}
	Point( const std::string& line ) {
		LoadLine(line);
	}
	bool LoadLine( const std::string& line ) {
		return ( sscanf(line.c_str(), "%f %f", &lon, &lat) == 2 );
	}

   inline const T& Lat() const { return lat; }
   inline	T& Lat() { return lat; }
   inline const T& Lon() const { return lon; }
   inline	T& Lon() { return lon; }
   friend std::ostream& operator << (std::ostream& o, Point a) { 
      //o << "(" << a.lon << ", " << a.lat<< ")"; 
      o.setf(std::ios::fixed);
      o << std::left << std::setprecision(4) << a.lon << " " << a.lat; 
      return o; 
   }
};
#endif


template < class T >
class DataPoint : public Point<T> {
   T dis, data;

public:
   DataPoint(const Point<T> ptin, const T datain = -12345., const T disin = -12345. )
      : Point<T>(ptin)
      , dis(disin), data(datain) {}

   DataPoint(const T lonin = -12345., const T latin = -12345., const T datain = -12345., const T disin = -12345. )
      : Point<T>(lonin, latin)
      , dis(disin), data(datain) {}

   inline const T& Dis() const { return dis; }
   inline	T& Dis() { return dis; }
   inline const T& Data() const { return data; }
   inline	T& Data() { return data; }

   friend std::ostream& operator << (std::ostream& o, DataPoint a) { 
      o << a.lon << " " << a.lat << " " <<a.data<<" "<<a.dis; 
      return o; 
   }

};


class Map {
public:
   //Map( const char *inname );
   Map( const char *inname, const float grdlon = 1.0, const float grdlat = 1.0 );
   Map( const char *inname, const Point<float>& srcin, const float grdlon = 1.0, const float grdlat = 1.0 );
   Map( const Map& );
   Map( Map&& );
   Map& operator= ( const Map& );
   Map& operator= ( Map&& );
   ~Map();

	/* ------------ set source location ------------ */
	//void SetSource( const Point<float>& src );

	/* ----- map boundaries ----- */
	float LonMin() const;
	float LonMax() const;
	float LatMin() const;
	float LatMax() const;

	/* ------------ compute number of points near the given location ------------ */
	float NumberOfPoints(Point<float> rec, const float xhdis, const float yhdis) const {
		float loneff, lateff;
		NumberOfPoints( rec, xhdis, yhdis, loneff, lateff );
	}
	float NumberOfPoints(Point<float> rec, const float xhdis, const float yhdis, float& loneff, float& lateff) const;

   /* ------------ compute average value on the point rec ------------ */
   float PointAverage(Point<float> Prec, float hdis ) {
      float weit;
      return PointAverage(Prec, hdis, weit);
   }
   float PointAverage(Point<float> rec, float hdis, float& weit);

   /* ------------ compute average value along the path src-rec ------------ */
   float PathAverage(Point<float> Prec, float lamda) {
      float perc;
      return (PathAverage(Prec, lamda, perc)).Data();
   }
   DataPoint<float> PathAverage(Point<float> Prec, float lamda, float& perc);

   /* ------------ compute average along the path src-rec weighted by the reciprocal of map values ------------ */
   float PathAverage_Reci(Point<float> Prec, float lamda) {
      float perc;
      return (PathAverage_Reci(Prec, lamda, perc)).Data();
   }
   DataPoint<float> PathAverage_Reci(Point<float> Prec, float lamda, float& perc, const std::string outname = "");
  
protected:
	static constexpr float NaN = -12345.;

private:
   struct Mimpl;
   std::unique_ptr<Mimpl> pimplM;
/*
   std::string fname;
   Point<float> src;
   float dismax, disamax;
   std::vector< DataPoint<float> > dataV;
   Array2D< DataPoint<float> > dataM;
*/
};


#endif
