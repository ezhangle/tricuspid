/******************************************************/
/*                                                    */
/* tricuspid.cpp - a curve with three cusps per side  */
/*                                                    */
/******************************************************/
/* Copyright 2018,2019 Pierre Abbat.
 * This file is part of Tricuspid.
 * 
 * Tricuspid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Tricuspid is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Tricuspid. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <bezitopo.h>
using namespace std;

void testcircle()
{
  Circle unit(xy(0,0),1.);
  Circle triple(xy(0,0),3.);
  Circle xaxis(xy(0,0),0,0);
  Circle yaxis(xy(0,0),DEG90,0);
  PostScript ps;
  int i,shortCount=0,longCount=0;
  xyz sta1,sta3;
  vector<segment> lines;
  polyline envelope;
  ps.open("circle.ps");
  ps.setpaper(papersizes["A4 portrait"],0);
  ps.prolog();
  ps.startpage();
  ps.setscale(-3,-3,3,3,degtobin(0));
  for (i=0;i<1080;i+=5)
  {
    sta1=unit.station(degtorad(i));
    sta3=triple.station(degtorad(i));
    ps.line2p(sta1,sta3);
    if (dist(sta1,sta3)<3.16227733)
      shortCount++;
    if (dist(sta1,sta3)>3.162278) // The distance at 45°/135° is sqrt(10).
      longCount++;
    lines.push_back(segment(sta1,sta3));
  }
  ps.spline(unit.approx3d(0.1/ps.getscale()));
  ps.spline(triple.approx3d(0.1/ps.getscale()));
  ps.setcolor(1,0,1);
  for (i=0;i<lines.size();i++)
    envelope.insert(intersection(lines[i],lines[(i+1)%lines.size()]));
  ps.spline(envelope.approx3d(0.1/ps.getscale()));
  ps.endpage();
  ps.startpage();
  ps.setscale(-1,-1,1,1,degtobin(0));
  ps.setcolor(0,0,0);
  for (i=8388608-DEG180;i<=DEG180;i+=16777216)
    ps.line2p(xaxis.station(tan(i)),yaxis.station(cot(i)));
  ps.spline(xaxis.approx3d(0.1/ps.getscale()));
  ps.spline(yaxis.approx3d(0.1/ps.getscale()));
  ps.endpage();
  ps.close();
}

xy curvePoint(int angle,int a,int b)
/* Computes the point on the curve defined by circles of sizes 1/a and 1/b,
 * going a distance bintorad(angle) around both circles.
 * a, b, a+b, and a-b must all be nonzero:
 * a==0 or b==0 means the circle is infinite.
 * a-b==0 means that the line goes through the same point twice and is indeterminate.
 * a+b==0 means that all the lines are parallel, so the envelope is at infinity.
 */
{
  xy apoint,bpoint,aspeed,bspeed;
  double atorque,btorque;
  apoint=cossin(angle*a)/a;
  bpoint=cossin(angle*b)/b;
  aspeed=cossin(angle*a+DEG90);
  bspeed=cossin(angle*b+DEG90);
  atorque=dot(turn90(aspeed),apoint-bpoint);
  btorque=dot(turn90(bspeed),bpoint-apoint);
  return (apoint*btorque+bpoint*atorque)/(atorque+btorque);
}

xy curvePoint(double angle,double a,double b)
/* Same as above, but the b circle is rotated 180° so that a cusp, not an
 * asymptote, is at angle=0. angle is in radians.
 */
{
  xy apoint,bpoint,aspeed,bspeed;
  double atorque,btorque;
  apoint=cossin(angle*a)/a;
  bpoint=-cossin(angle*b)/b;
  aspeed=cossin(angle*a+M_PI/2);
  bspeed=-cossin(angle*b+M_PI/2);
  atorque=dot(turn90(aspeed),apoint-bpoint);
  btorque=dot(turn90(bspeed),bpoint-apoint);
  return (apoint*btorque+bpoint*atorque)/(atorque+btorque);
}

polyline curvePart(int startAngle,int endAngle,int a,int b)
{
  int i,n;
  double h;
  xy pnt;
  polyline ret;
  const int clipto=256;
  int midAngle=startAngle+(endAngle-startAngle)/2;
  brent br;
  startAngle=br.init(startAngle-midAngle,curvePoint(startAngle,a,b).length()-clipto,
		     0,curvePoint(midAngle,a,b).length()-clipto,true)+midAngle;
  while (!br.finished())
    startAngle=br.step(curvePoint(startAngle,a,b).length()-clipto)+midAngle;
  endAngle=br.init(endAngle-midAngle,curvePoint(endAngle,a,b).length()-clipto,
		   0,curvePoint(midAngle,a,b).length()-clipto,true)+midAngle;
  while (!br.finished())
    endAngle=br.step(curvePoint(endAngle,a,b).length()-clipto)+midAngle;
  h=(endAngle-startAngle)/rint((endAngle-startAngle)/DEG1);
  n=lrint((endAngle-startAngle)/h);
  for (i=0;i<=n;i++)
  {
    pnt=curvePoint(startAngle+lrint(i*h),a,b);
    ret.insert(pnt);
  }
  ret.open();
  return ret;
}

polyline curvePart(double startAngle,double endAngle,double a,double b)
{
  int i,n;
  double h;
  xy pnt;
  polyline ret;
  const int clipto=256;
  double midAngle=(startAngle+endAngle)/2;
  brent br;
  startAngle=br.init(startAngle-midAngle,curvePoint(startAngle,a,b).length()-clipto,
		     0,curvePoint(midAngle,a,b).length()-clipto,false)+midAngle;
  while (!br.finished())
    startAngle=br.step(curvePoint(startAngle,a,b).length()-clipto)+midAngle;
  endAngle=br.init(endAngle-midAngle,curvePoint(endAngle,a,b).length()-clipto,
		   0,curvePoint(midAngle,a,b).length()-clipto,false)+midAngle;
  while (!br.finished())
    endAngle=br.step(curvePoint(endAngle,a,b).length()-clipto)+midAngle;
  h=(endAngle-startAngle)/rint((endAngle-startAngle)*57);
  n=lrint((endAngle-startAngle)/h);
  for (i=0;i<=n;i++)
  {
    pnt=curvePoint(startAngle+i*h,a,b);
    ret.insert(pnt);
  }
  ret.open();
  return ret;
}

void drawcurve(int a,int b,PostScript &ps)
{
  Circle aCircle(xy(0,0),1./a);
  Circle bCircle(xy(0,0),1./b);
  int nparts=abs(a-b);
  int i;
  polyline part;
  ps.setpaper(papersizes["A4 portrait"],0);
  ps.prolog();
  ps.startpage();
  ps.setscale(-1,-1,1,1,degtobin(0));
  ps.setcolor(1,0,1);
  ps.spline(aCircle.approx3d(0.1/ps.getscale()));
  ps.spline(bCircle.approx3d(0.1/ps.getscale()));
  ps.setcolor(0,0,0);
  for (i=0;i<nparts;i++)
  {
    part=curvePart(radtobin(i*2*M_PI/nparts)+FURMAN1,radtobin((i+1)*2*M_PI/nparts)-FURMAN1,a,b);
    ps.spline(part.approx3d(0.1/ps.getscale()));
  }
  ps.endpage();
}

void drawcurve(double a,double b,PostScript &ps)
{
  int i,n;
  double h=degtorad(10);
  Circle aCircle(xy(0,0),1./a);
  Circle bCircle(xy(0,0),1./b);
  double nparts=abs(a-b);
  polyline part;
  ps.setpaper(papersizes["A4 portrait"],0);
  ps.prolog();
  ps.startpage();
  ps.setscale(-1,-1,1,1,degtobin(0));
  ps.setcolor(1,0,1);
  ps.spline(aCircle.approx3d(0.1/ps.getscale()));
  ps.spline(bCircle.approx3d(0.1/ps.getscale()));
  n=trunc(M_PI/nparts/h);
  ps.setcolor(0,1,0);
  for (i=-n;i<=n;i++)
    ps.line2p(-aCircle.station(i*h),bCircle.station(i*h));
  ps.setcolor(0,0,0);
  // Unlike the integer case, it makes no sense to draw multiple parts if a and b are irrational.
  part=curvePart(-M_PI/nparts+1e-9,M_PI/nparts-1e-9,a,b);
  ps.spline(part.approx3d(0.1/ps.getscale()));
  ps.endpage();
}

bool check142(int angle)
/* The curves for (-1,2) and (1,4) look identical, except that one is -2 times
 * as big as the other. Check that the points on circles of radius 1, 1/4,
 * and -1/2 are in a straight line.
 */
{
  xy pt1=cossin(angle);
  xy pt4=cossin(angle*4)/4;
  xy pt2=cossin(-angle*2)/-2;
  return fabs(pldist(pt4,pt2,pt1))<1e-6;
}

double est2Deriv(double a,double b)
/* Estimates the second derivative of x with respect to the angle at the cusp.
 * The cusp changes from triple to single when the second derivative is 0.
 * When it is 0, the cusp is a 4/5 cusp, instead of the usual 2/3 cusp(s).
 */
{
  const int num=26;
  double xval[num],x0,der2[num],angle;
  int i;
  x0=curvePoint(0.,a,b).getx();
  for (i=0,angle=1;i<num;i++,angle/=2)
  {
    xval[i]=curvePoint(angle,a,b).getx();
    der2[i]=2*(xval[i]-x0)/sqr(angle);
    //cout<<setw(2)<<i<<' '<<ldecimal(der2[i])<<endl;
  }
  return der2[14];
}

int main(int argc, char *argv[])
{
  int a,b;
  double bd;
  brent br;
  PostScript ps;
  testcircle();
  ps.open("tricuspid.ps");
  for (b=2;b<6;b++)
    for (a=-b;a<b;a++)
      if (gcd(abs(a),b)==1)
	drawcurve(a,b,ps);
  for (a=1,b=2;b<35;b+=a,a=b-a)
    drawcurve(b-a,b,ps);
  cout<<check142(DEG1)<<endl;
  drawcurve(1,3,ps); //in
  drawcurve(1,4,ps); //out
  drawcurve(2,7,ps); //in
  drawcurve(3,10,ps); //in
  drawcurve(1.,10/3.,ps);
  drawcurve(1.,3.5,ps);
  est2Deriv(1.,10/3);
  bd=br.init(3.3,est2Deriv(1.,3.3),4,est2Deriv(1.,4.));
  while (!br.finished())
    bd=br.step(est2Deriv(1.,bd));
  cout<<"4/5 cusp at b="<<bd<<endl;
  drawcurve(1.,bd,ps);
  est2Deriv(1.,bd);
  return 0;
}
