#ifndef DISAZI_H
#define DISAZI_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
//using namespace std;

#ifndef POINT
#define POINT
template <class T>
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
   inline       T& Lat() { return lat; }
   inline const T& Lon() const { return lon; }
   inline       T& Lon() { return lon; }
   friend std::ostream& operator << (std::ostream& o, Point a) {
      //o << "(" << a.lon << ", " << a.lat<< ")"; 
      o.setf(std::ios::fixed);
      o << std::left << std::setprecision(4) << a.lon << " " << a.lat;
      return o;
   }
};
#endif


template <class T>
class Path {
public:
   Path(const T lo1in=-12345., const T la1in=-12345., 
        const T lo2in=-12345., const T la2in=-12345.)
      : long1(lo1in), lati1(la1in), long2(lo2in), lati2(la2in)
      , dist(-12345.), alpha1(-12345.), alpha2(-12345.) {}

   Path(const Point<T>& p1, const Point<T>& p2)
      : long1(p1.Lon()), lati1(p1.Lat()), long2(p2.Lon()), lati2(p2.Lat())
      , dist(-12345.), alpha1(-12345.), alpha2(-12345.) {}

   ~Path(){}

	Point<T> P1() { return Point<T>(long1, lati1); }
	Point<T> P2() { return Point<T>(long2, lati2); }
	Point<T> PC() { return Point<T>(0.5*(long1+long2), 0.5*(lati1+lati2)); }

	friend std::ostream& operator << ( std::ostream& o, class Path p ) {
		o<<p.long1<<" "<<p.lati1<<"   "<<p.long2<<" "<<p.lati2;
		return o;
	}

   inline const T& Dist() {
      if( dist == -12345. ) calc_dist();
      return dist;
   }

   inline const T& Azi1() {
      if( alpha1 == -12345. ) calc_azimuth();
      return alpha1;
   }

   inline const T& Azi2() {
      if( alpha2 == -12345. ) calc_azimuth();
      return alpha2;
   }

private:
   T dist, alpha1, alpha2;
   T lati1, long1, lati2, long2;

   void calc_dist() {
      if( lati1==-12345. || long1==-12345. || lati2==-12345. || long2==-12345. )
	 throw std::runtime_error("Error(calc_azimuth): empty location(s)!");

      int i,flag=0;
      T latio1,longo1,dlt_long;//dlt_lati;
      T Rx;//R;
      T U1,U2;
      T Ds,alpha2;
      //T ctr_angl;
      T cv,cv1,cv2,cv3,cv4,cv5,cvC,numda1;
      T numda;
      T mius, cvA, cvB, deltacv;

      T pi = 3.1415926535898; //4.0*atan(1.0);
      T pio180 = 0.0174532925199433;
      T Ra = 6378.137, Rb = 6356.7523142;
      T f = 0.00335281066474748; //1/298.257223563;

   begin:
      //long1=long1-(int)floor(long1)/360*360;
      //if(long1<0.) long1+=360.;
      while( long1 < 0. ) long1 += 360.;
      while( long1 >= 360. ) long1 -= 360.;
      //long2=long2-(int)floor(long2)/360*360;
      //if(long2<0.) long2+=360.;
      while( long2 < 0. ) long2 += 360.;
      while( long2 >= 360. ) long2 -= 360.;
      if(lati1==-lati2 && fabs(long1-long2)==180){
	 dist = 20003.917357;
	 return;
      }
      //dlt_lati=fabs(lati2-lati1);
      dlt_long=long2-long1;

      if (dlt_long > 180.000)  dlt_long = 360.000000 - dlt_long;
      if (dlt_long < -180.000) dlt_long = 360.000 - fabs(dlt_long);
      dlt_long = fabs(dlt_long);

      U1 = atan((1-f)*tan(lati1*pio180));
      U2 = atan((1-f)*tan(lati2*pio180));
      T sinU1 = sin(U1), cosU1 = cos(U1);
      T sinU2 = sin(U2), cosU2 = cos(U2);
      dlt_long = dlt_long*pio180;

      numda = dlt_long;
      numda1 = numda;
      i=0;
      do {
	 i++;
	 numda = numda1;
	 T sinnumda = sin(numda), cosnumda = cos(numda);
	 T ftmp1 = cosU2*sinnumda, ftmp2 = cosU1*sinU2-sinU1*cosU2*cosnumda;
	 cv1 =  sqrt( ftmp1*ftmp1 + ftmp2*ftmp2 );
	 cv2 = sinU1*sinU2+ cosU1*cosU2*cosnumda;
	 cv = atan2(cv1,cv2);
	 if(cv==0) cv = 1.e-10;//0.00000000001;
	 cv3 = cosU1*cosU2*sinnumda/sin(cv);
	 cv4 = 1 - cv3*cv3;
	 if(cv4==0) cv4 = 1.e-10;//0.00000000001;
	 cv5 = cos(cv) - 2*sinU1*sinU2/cv4;
	 cvC = f*0.0625*cv4*(4 + f*(4 - 3*cv4));
	 numda1 = dlt_long + (1-cvC)*f*cv3*(cv + cvC*cv1*(cv5 + cvC*cv2*(-1 +2*cv5*cv5)));
	 if(i>50){
	    flag=1;
	    latio1=lati1*pio180;
	    ftmp1 = Ra*cos(latio1); ftmp2 = Rb*sin(latio1);
	    T ftmp3 = Ra*ftmp1, ftmp4 = Rb*ftmp2;
	    Rx=sqrt((ftmp3*ftmp3+ftmp4*ftmp4)/(ftmp1*ftmp1+ftmp2*ftmp2));
	    Rx=pi*(Rx+Ra)*0.5;
	    latio1=lati1;
	    longo1=long1;
	    lati1=lati2;
	    long1=long2;
	    lati2=-latio1;
	    long2=longo1+180;
	    goto begin;
	 }
      } while (fabs(numda - numda1) > 1e-5);
      mius = cv4*(Ra*Ra/(Rb*Rb) - 1.);
      cvA = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)));
      cvB = mius/1024.*(256+ mius*(-128 + mius*(74 - 47*mius)));
      deltacv = cvB*cv1*(cv5 +cvB/4*(cv2*(-1 + 2*cv5*cv5)-cvB/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ));
      dist = Rb * cvA *(cv - deltacv);

      if(flag==1){
	 alpha1=atan2(cosU2*sin(dlt_long), cosU1*sinU2 - sinU1*cosU2*cos(dlt_long))/pio180;
	 alpha2=atan2(cosU1*sin(dlt_long), -sinU1*cosU2 + cosU1*sinU2*cos(dlt_long))/pio180;
	 if( fabs(long2-long1)>180 ) { alpha1 = 360.-alpha1; alpha2 = 360.-alpha2;}
	 if( long2 < long1 ) { alpha1 = 360.-alpha1; alpha2 = 360.-alpha2; }
	 float theta = alpha1;
	 if(theta>180) theta=360-theta;
	 theta=fabs(90-theta);
	 Ds=Rx*(90-theta)/90.+(Ra+Rb)*pio180*theta;
	 dist=Ds-dist;
      }

   }

   void calc_azimuth() {
      if( lati1==-12345. || long1==-12345. || lati2==-12345. || long2==-12345. )
	 throw std::runtime_error("Error(calc_azimuth): empty location(s)!");

      T dlt_long;//dlt_lati;
      T U1,U2;

      T pi = 3.14159265358979324; //4.0*atan(1.0);
      T pio180 = 0.017453292519943295;
      //Ra = 6378.137;
      //Rb = 6356.7523142;
      T f = 0.0033528106647474805; //1/298.257223563;
      long1=long1-(int)floor(long1)/360*360;
      if(long1<0) long1+=360;
      long2=long2-(int)floor(long2)/360*360;
      if(long2<0) long2+=360;
      if(lati1==-lati2 && fabs(long1-long2)==180){
         alpha1 = alpha2 = 999.;
	 return;
      }
      //dlt_lati=fabs(lati2-lati1);
      dlt_long=long2-long1;

      if (dlt_long > 180.000)  dlt_long = 360.000000 - dlt_long;
      if (dlt_long < -180.000) dlt_long = 360.000 - fabs(dlt_long);
      dlt_long = fabs(dlt_long);

      U1 = atan((1-f)*tan(lati1*pio180));
      U2 = atan((1-f)*tan(lati2*pio180));
      T sinU1 = sin(U1), sinU2 = sin(U2);
      T cosU1 = cos(U1), cosU2 = cos(U2);
      dlt_long = dlt_long*pi/180;

      alpha1=atan2(cosU2*sin(dlt_long), cosU1*sinU2 - sinU1*cosU2*cos(dlt_long))*180/pi;
      alpha2=atan2(cosU1*sin(dlt_long), -sinU1*cosU2 + cosU1*sinU2*cos(dlt_long))*180/pi;
      if( fabs(long2-long1)>180 ) { alpha1 = 360 - alpha1; alpha2 = 360 - alpha2; }
      if( long2 < long1 ) { alpha1 = 360 - alpha1; alpha2 = 360 - alpha2; }

   }

};

/*
template <class T>
bool Check_azi_cov(T *lon, T *lat, int nsta, T mdis, T *flag)
{
   int i, j, jj, k, ndata;
   T azi[300], azimin, azitmp, dist;
   T mazi=110.;

   for(i=0;i<nsta;i++){
      if(flag[i]==-1) continue;
      ndata=0;
      for(j=0;j<nsta;j++){
         if(j==i || flag[j]==-1) continue;
         calc_dist(lat[i], lon[i], lat[j], lon[j], &dist);
         if(dist>mdis) continue;
         calc_azimuth(lat[i], lon[i], lat[j], lon[j], &azi[ndata]);
         ndata++;
      }
      //cout<<"Check azi cov: "<<ndata<<" nearby stations"<<endl;
      for(j=0;j<ndata;j++) {
         azimin=azi[j]; jj=j;
         for(k=j+1;k<ndata;k++)
            if(azimin>azi[k]) {
              azimin=azi[k];
              jj=k;
            }
         if(jj==j) continue;
         azitmp=azi[j];
         azi[j]=azi[jj];
         azi[jj]=azitmp;
      }
      if(azi[0]+360.-azi[ndata-1]>mazi) {
         flag[i]=-2;
         continue;
      }
      for(j=1;j<ndata;j++)
         if(azi[j]-azi[j-1]>mazi) {
            flag[i]=-2;
            break;
         }
   }
   for(i=0;i<nsta;i++) if(flag[i]==-2) flag[i]=-1;
   return 1;
}

*/

#endif
