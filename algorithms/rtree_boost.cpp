#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

// to store queries results
#include <vector>
#include <time.h>
#include <stdlib.h>

// just for output
#include <iostream>
#include <boost/foreach.hpp>

using namespace std;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

FILE* fin;
FILE* fout;

int main()
{
    fin=fopen("dataobs.txt","r");
    if(fin==NULL)
		printf("Error opening file");

    typedef bg::model::point<double, 2, bg::cs::cartesian> point;
    typedef bg::model::box<point> box;
    typedef std::pair<box, unsigned> value;

    // create the rtree using default constructor
    bgi::rtree< value, bgi::quadratic<100001> > rtree;

    // create some values
    for ( unsigned i = 0 ; i < 100000 ; ++i )
    {
        double xl,xh,yl,yh;
        fscanf(fin," %lf %lf %lf %lf",&xl,&xh,&yl,&yh);
        // create a box
        box b(point(xl, yl), point(xh, yh));
        //cout<<i<<endl;
        // insert new value
        rtree.insert(std::make_pair(b, i));
    }

    fclose(fin);



    fout=fopen("rq.txt","w+");
    if(fout==NULL)
		printf("Error opening file");


    fin=fopen("datadps.txt","r");

	if(fin==NULL)
		printf("Error opening file");

    const clock_t t1 = clock();
    srand(time(NULL));

    // find values intersecting some area defined by a box
    int i,j;
    for(i=0;i<256;i++)
    {
        double px,py;
        fscanf(fin," %lf %lf",&px,&py);
        fprintf(fout,"%lf %lf\n",px,py);
        box query_box(point(px-40.0,py-40.0), point(px+40.0,py+40.0));
        std::vector<value> result_s;
        rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));
        fprintf(fout,"%d\n",result_s.size());
        for(j=0;j<result_s.size();j++)
        {
            point minn=result_s[j].first.min_corner();
            point maxx=result_s[j].first.max_corner();

            double xl=minn.get<0>();
            double xh=maxx.get<0>();
            double yl=minn.get<1>();
            double yh=maxx.get<1>();
            fprintf(fout,"%lf %lf %lf %lf\n",xl,xh,yl,yh);
        }
        fprintf(fout,"\n");

    }

    fclose(fin);
    fclose(fout);





    const clock_t t2 = clock();
    cout<<t1<<" "<<t2<<" "<<t2-t1<<endl;

    // note: in Boost.Geometry WKT representation of a box is polygon

    // display results


//    BOOST_FOREACH(value const& v, result_s)
//    {
//        point minn=v.first.min_corner();
//        point maxx=v.first.max_corner();
//
//        double xl=minn.get<0>();
//        double xh=maxx.get<0>();
//        double yl=minn.get<1>();
//        double yh=maxx.get<1>();
//
//        //std::cout <<xl<<"\t"<<xh<<"\t"<<yl<<"\t"<<yh<< std::endl;
//
//    }
        //std::cout << bg::wkt<box>(v.first) << std::endl;
        //std::cout << bg::wkt<box>(v.first) << " - " << v.second << std::endl;



    return 0;
}
