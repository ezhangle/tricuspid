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

int main(int argc, char *argv[])
{
  cout<<"tricuspid\n";
  testcircle();
  return 0;
}
