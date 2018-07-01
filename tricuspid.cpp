/******************************************************/
/*                                                    */
/* tricuspid.cpp - a curve with three cusps per side  */
/*                                                    */
/******************************************************/

#include <iostream>
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

polyline curvePart(int startAngle,int endAngle,int a,int b)
{
  int i,n;
  double h;
  xy pnt;
  polyline ret;
  h=(endAngle-startAngle)/rint((endAngle-startAngle)/DEG1);
  n=lrint((endAngle-startAngle)/h);
  for (i=0;i<=n;i++)
  {
    pnt=curvePoint(startAngle+lrint(i*h),a,b);
    if (pnt.length()<256)
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

int main(int argc, char *argv[])
{
  int a,b;
  PostScript ps;
  testcircle();
  ps.open("tricuspid.ps");
  for (b=2;b<6;b++)
    for (a=-b;a<b;a++)
      if (gcd(abs(a),b)==1)
	drawcurve(a,b,ps);
  for (a=1,b=2;b<2;b+=a,a=b-a)
    drawcurve(b-a,b,ps);
  return 0;
}
