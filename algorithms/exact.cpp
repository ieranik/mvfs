#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <map>

#include <windows.h>
#include <glut.h>

#define pi (2*acos(0.0))
#define eps 0.0000001

using namespace std;

double max2(double a,double b)
{
    return (a>b)?a:b;
}

double min2(double a,double b)
{
    return (a<b)?a:b;
}

class point
{
public:
    double x,y;
    point(double xx=0.0, double yy=0.0)
    {
        x=xx;
        y=yy;
    }
};

class triangle
{
public:
    point a, b, c;
    unsigned int bitmap;
    int cnt;
    double area;
    triangle(point aa=point(), point bb=point(), point cc=point(),unsigned int bm=0,int cn=0)
    {
        a=aa;
        b=bb;
        c=cc;
        bitmap=bm;
        cnt=cn;
    }
};

class triedge
{
public:
    point a, b;
    triedge(point aa=point(), point bb=point())
    {
        a=aa;
        b=bb;
    }
};

double tempcx, tempcy;

bool cmptrip(point a,point b)
{
    double aa=atan2(a.y-tempcy,a.x-tempcx);
    double ab=atan2(b.y-tempcy,b.x-tempcx);
    return (aa<ab);
}

bool samepoint(point a, point b)
{
    if(fabs(a.x-b.x)<eps&&fabs(a.y-b.y)<eps)return true;
    return false;
}

double areat(triangle t)
{
    return t.a.x*(t.b.y-t.c.y)+t.b.x*(t.c.y-t.a.y)+t.c.x*(t.a.y-t.b.y);
}

double area3p(point a, point b, point c)
{
    return a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y);
}

bool pinsidet(point p, triangle t)
{
    if(area3p(t.a,t.b,t.c)*area3p(t.a,t.b,p)>-eps&&area3p(t.a,t.b,t.c)*area3p(t.a,p,t.c)>-eps&&area3p(t.a,t.b,t.c)*area3p(p,t.b,t.c)>-eps)
        return true;
    return false;
}

bool psinsidet(point p, triangle t)
{
    if(area3p(t.a,t.b,t.c)*area3p(t.a,t.b,p)>0&&area3p(t.a,t.b,t.c)*area3p(t.a,p,t.c)>0&&area3p(t.a,t.b,t.c)*area3p(p,t.b,t.c)>0)
        return true;
    return false;
}

bool tinsidet(point a,point b,point c,triangle t)
{
    if(pinsidet(a,t)&&pinsidet(b,t)&&pinsidet(c,t))return true;
    return false;
}

bool doesint(point a, point b, point c, point d)
{
    if(area3p(a,b,c)*area3p(a,b,d)<=0&&area3p(c,d,a)*area3p(c,d,b)<=0)
        return true;
    return false;
}

point intp(point a, point b, point c, point d)
{
    double cc=fabs(area3p(a,b,c));
    double dd=fabs(area3p(a,b,d));
    double r=cc/(cc+dd);
    double x=c.x*(1-r)+d.x*r;
    double y=c.y*(1-r)+d.y*r;
    return point(x,y);
}

bool ifaonbc(point a, point b, point c)
{
    if(fabs((c.y-a.y)*(a.x-b.x)-(a.y-b.y)*(c.x-a.x))<eps)return true;
    return false;
}

bool doesinthard(point a, point b, point c, point d)
{
    if(area3p(a,b,c)*area3p(a,b,d)<-eps&&area3p(c,d,a)*area3p(c,d,b)<-eps)
        return true;
    return false;
}


bool hardint(point a,point b,point c,vector<triedge> segs)
{
    int i;
    for(i=0;i<segs.size();i++)
    {
        if(doesinthard(a,b,segs[i].a,segs[i].b))return true;
        if(doesinthard(b,c,segs[i].a,segs[i].b))return true;
        if(doesinthard(c,a,segs[i].a,segs[i].b))return true;
    }
    return false;
}


void pp(point p)
{
    printf("%lf %lf\n",p.x,p.y);
}

void pt(triangle t)
{
    pp(t.a);
    pp(t.b);
    pp(t.c);
    cout<<endl;
}

point cog(triangle t)
{
    point ret;
    ret.x=(t.a.x+t.b.x+t.c.x)/3.0;
    ret.y=(t.a.y+t.b.y+t.c.y)/3.0;
    return ret;
}

vector<triangle> removedup(vector<triangle> vt)
{
    vector<triangle> ret;
    vector<point> c;
    int i,j;
    //cout<<vt.size()<<endl<<endl;
    //for(i=0;i<vt.size();i++)pt(vt[i]);
    for(i=0;i<vt.size();i++)
    {
        bool isfound=false;
        point cg=cog(vt[i]);
        for(j=0;j<c.size();j++)
        {
            if(samepoint(cg,c[j]))
            {
                isfound=true;
                break;
            }
        }
        if(isfound==false)
        {
            c.push_back(cg);
            //if(fabs(areat(vt[i]))>eps)
                ret.push_back(vt[i]);
        }
    }
    //cout<<ret.size();
    return ret;
}

vector< vector<triangle> > inttri(triangle a,triangle b)
{
    vector< vector<triangle> > ret;
    vector<triangle> aandb,anotb,bnota;

    vector<point> plista,plistb;
    vector<triedge> segs;

    int i,j;


    while(true)
    {
        if(fabs(area3p(a.a,a.b,b.a))<eps||fabs(area3p(a.a,a.b,b.b))<eps||fabs(area3p(a.a,a.b,b.c))<eps||fabs(area3p(a.b,a.c,b.a))<eps||fabs(area3p(a.b,a.c,b.b))<eps||fabs(area3p(a.b,a.c,b.c))<eps||fabs(area3p(a.c,a.a,b.a))<eps||fabs(area3p(a.c,a.a,b.b))<eps||fabs(area3p(a.c,a.a,b.c))<eps)
        {
            int ri;
            double rd;
            srand(time(NULL));
            ri=rand()%200-200;
            rd=((double)ri)*eps;
            a.a.x+=rd;
            ri=rand()%200-200;
            rd=((double)ri)*eps;
            a.a.y+=rd;
            ri=rand()%200-200;
            rd=((double)ri)*eps;
            a.b.x+=rd;
            ri=rand()%200-200;
            rd=((double)ri)*eps;
            a.b.y+=rd;
            ri=rand()%200-200;
            rd=((double)ri)*eps;
            a.c.x+=rd;
            ri=rand()%200-200;
            rd=((double)ri)*eps;
            a.c.y+=rd;
        }
        else break;
    }

    bool didint=false;

    if(doesint(a.a,a.b,b.a,b.b)){
        didint=true;
        plista.push_back(intp(a.a,a.b,b.a,b.b));
        plistb.push_back(intp(a.a,a.b,b.a,b.b));
    }

    if(doesint(a.a,a.b,b.b,b.c)){
        didint=true;
        plista.push_back(intp(a.a,a.b,b.b,b.c));
        plistb.push_back(intp(a.a,a.b,b.b,b.c));
    }

    if(doesint(a.a,a.b,b.c,b.a)){
        didint=true;
        plista.push_back(intp(a.a,a.b,b.c,b.a));
        plistb.push_back(intp(a.a,a.b,b.c,b.a));
    }



    if(doesint(a.b,a.c,b.a,b.b)){
        didint=true;
        plista.push_back(intp(a.b,a.c,b.a,b.b));
        plistb.push_back(intp(a.b,a.c,b.a,b.b));
    }

    if(doesint(a.b,a.c,b.b,b.c)){
        didint=true;
        plista.push_back(intp(a.b,a.c,b.b,b.c));
        plistb.push_back(intp(a.b,a.c,b.b,b.c));
    }

    if(doesint(a.b,a.c,b.c,b.a)){
        didint=true;
        plista.push_back(intp(a.b,a.c,b.c,b.a));
        plistb.push_back(intp(a.b,a.c,b.c,b.a));
    }



    if(doesint(a.c,a.a,b.a,b.b)){
        didint=true;
        plista.push_back(intp(a.c,a.a,b.a,b.b));
        plistb.push_back(intp(a.c,a.a,b.a,b.b));
    }

    if(doesint(a.c,a.a,b.b,b.c)){
        didint=true;
        plista.push_back(intp(a.c,a.a,b.b,b.c));
        plistb.push_back(intp(a.c,a.a,b.b,b.c));
    }

    if(doesint(a.c,a.a,b.c,b.a)){
        didint=true;
        plista.push_back(intp(a.c,a.a,b.c,b.a));
        plistb.push_back(intp(a.c,a.a,b.c,b.a));
    }

    if(didint==false)
    {
        if(pinsidet(a.a,b)==false && pinsidet(b.a,a)==false)
        {
            anotb.push_back(a);
            bnota.push_back(b);
            ret.push_back(aandb);
            ret.push_back(anotb);
            ret.push_back(bnota);
            //printf("%d %d %d\n",aandb.size(),anotb.size(),bnota.size());
            return ret;
        }
        if(pinsidet(a.a,b)==true)
        {
            aandb.push_back(a);
        }
        else aandb.push_back(b);
    }

    plista.push_back(a.a);
    //plistb.push_back(a.a);

    plista.push_back(a.b);
    //plistb.push_back(a.b);

    plista.push_back(a.c);
    //plistb.push_back(a.c);


    //plista.push_back(b.a);
    plistb.push_back(b.a);

    //plista.push_back(b.b);
    plistb.push_back(b.b);

    //plista.push_back(b.c);
    plistb.push_back(b.c);


    tempcx=(a.a.x+a.b.x+a.c.x)/3.0;
    tempcy=(a.a.y+a.b.y+a.c.y)/3.0;

    sort(plista.begin(),plista.end(),cmptrip);

    //for(i=0;i<plista.size();i++)pp(plista[i]);
    //printf("\n");


    tempcx=(b.a.x+b.b.x+b.c.x)/3.0;
    tempcy=(b.a.y+b.b.y+b.c.y)/3.0;

    sort(plistb.begin(),plistb.end(),cmptrip);

    //for(i=0;i<plistb.size();i++)pp(plistb[i]);


    segs.push_back(triedge(a.a,a.b));
    segs.push_back(triedge(a.b,a.c));
    segs.push_back(triedge(a.c,a.a));

    segs.push_back(triedge(b.a,b.b));
    segs.push_back(triedge(b.b,b.c));
    segs.push_back(triedge(b.c,b.a));


    plista.push_back(plista[0]);
    plistb.push_back(plistb[0]);
    for(i=0;i<plista.size()-1;i++)
    {
        for(j=0;j<plistb.size()-1;j++)
        {
            point ta=plista[i];
            point tb=plista[i+1];
            point tc=plistb[j];
            if(fabs(area3p(ta,tb,tc))<0.1)continue;
            if(hardint(ta,tb,tc,segs))continue;
            bool ina=tinsidet(ta,tb,tc,a);
            bool inb=tinsidet(ta,tb,tc,b);

            if(ina==false&&inb==false)continue;
            triangle tt=triangle(ta,tb,tc);
            if(ina==true&&inb==true)aandb.push_back(tt);
            else if(ina==true)
            {
               if(psinsidet(b.a,tt)==false&&psinsidet(b.b,tt)==false&&psinsidet(b.c,tt)==false)
                    anotb.push_back(triangle(ta,tb,tc));
            }
            else
            {
                if(psinsidet(a.a,tt)==false&&psinsidet(a.b,tt)==false&&psinsidet(a.c,tt)==false)
                    bnota.push_back(triangle(ta,tb,tc));
            }
            segs.push_back(triedge(ta,tb));
            segs.push_back(triedge(tb,tc));
            segs.push_back(triedge(tc,ta));
        }
    }
    for(i=0;i<plistb.size()-1;i++)
    {
        for(j=0;j<plista.size()-1;j++)
        {
            point ta=plistb[i];
            point tb=plistb[i+1];
            point tc=plista[j];
            if(fabs(area3p(ta,tb,tc))<0.1)continue;
            if(hardint(ta,tb,tc,segs))continue;
            bool ina=tinsidet(ta,tb,tc,a);
            bool inb=tinsidet(ta,tb,tc,b);

            if(ina==false&&inb==false)continue;
            triangle tt=triangle(ta,tb,tc);
            if(ina==true&&inb==true)aandb.push_back(tt);
            else if(ina==true)
            {
                if(psinsidet(b.a,tt)==false&&psinsidet(b.b,tt)==false&&psinsidet(b.c,tt)==false)
                    anotb.push_back(triangle(ta,tb,tc));
            }
            else
            {
                if(psinsidet(a.a,tt)==false&&psinsidet(a.b,tt)==false&&psinsidet(a.c,tt)==false)
                    bnota.push_back(triangle(ta,tb,tc));
            }
            segs.push_back(triedge(ta,tb));
            segs.push_back(triedge(tb,tc));
            segs.push_back(triedge(tc,ta));
        }
    }

    aandb=removedup(aandb);
    anotb=removedup(anotb);
    bnota=removedup(bnota);

    ret.push_back(aandb);
    ret.push_back(anotb);
    ret.push_back(bnota);
    //for(i=0;i<aandb.size();i++)
    {
        //cout<<areat(aandb[i])<<endl;
        //pt(aandb[i]);
    }
    //printf("%d %d %d\n",aandb.size(),anotb.size(),bnota.size());
    return ret;

}


class record
{
public:
    unsigned int bitmap;
    double area;
    int cnt;
    record(unsigned int bm=0,double d=0,int ct=0)
    {
        bitmap=bm;
        area=d;
        cnt=ct;
    }
};

bool comrec(record a,record b)
{
    return (a.cnt<b.cnt);
}

int numone(unsigned int n)
{
    int cnt=0;
    while(n)
    {
        if(n%2==1)cnt++;
        n/=2;
    }
    return cnt;
}

class seg
{
public:
    int id;
    point p1, p2;
    bool hor,dir;
    seg(point p11=point(), point p22=point(),bool horr=false,bool direction=false,int idd=-1)
    {
        p1=p11;
        p2=p22;
        hor=horr;
        dir=direction;
        id=idd;
    }
};

class obs
{
public:
    double xl,xh,yl,yh;
    obs(double xxl=0.0, double xxh=0.0, double yyl=0.0, double yyh=0.0)
    {
        xl=xxl;
        xh=xxh;
        yl=yyl;
        yh=yyh;
    }
};

bool cmpxl(obs a,obs b)
{
	if(a.xl<b.xl)return true;
	return false;
}

double dist(point p1,point p2)
{
    return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}

bool cmpseg(seg a,seg b)
{
    double aa=atan2(a.p1.y,a.p1.x);
    double ab=atan2(b.p1.y,b.p1.x);
    if(fabs(aa-ab)<eps)
    {
        if(dist(point(),a.p1)<dist(point(),b.p1))return true;
        else return false;
    }
    else if(aa<ab)return true;
    else return false;
}

obs obsrot(obs o,double angle)
{
    double xl=cos(angle)*o.xl-sin(angle)*o.yl;
    double yl=sin(angle)*o.xl+cos(angle)*o.yl;
    double xh=cos(angle)*o.xh-sin(angle)*o.yh;
    double yh=sin(angle)*o.xh+cos(angle)*o.yh;
    return obs(min2(xl,xh),max2(xl,xh),min2(yl,yh),max2(yl,yh));
}

seg segrot(seg s,double angle)
{
    seg ret;
    ret.hor=s.hor;
    ret.p1.x=cos(angle)*s.p1.x-sin(angle)*s.p1.y;
    ret.p1.y=sin(angle)*s.p1.x+cos(angle)*s.p1.y;
    ret.p2.x=cos(angle)*s.p2.x-sin(angle)*s.p2.y;
    ret.p2.y=sin(angle)*s.p2.x+cos(angle)*s.p2.y;
    return ret;
}

double sangle(seg a)
{
    if(a.dir)return atan2(a.p1.y,a.p1.x);
    else return atan2(a.p2.y,a.p2.x);
}

double eangle(seg a)
{
    if(a.dir)return atan2(a.p2.y,a.p2.x);
    else return atan2(a.p1.y,a.p1.x);
}

double angle1(seg a)
{
    return atan2(a.p1.y,a.p1.x);
}

point intersect(seg s,double angle)
{
    double t=tan(angle);
    if(s.hor==true)return point(s.p1.y/t,s.p1.y);
    return point(s.p1.x,s.p1.x*t);
}

seg nearinter(vector<seg> s,double angle)
{
    seg n=s[0];
    double neardis=dist(point(),intersect(s[0],angle));
    int i;
    for(i=1;i<s.size();i++)
    {
        double tmp=dist(point(),intersect(s[i],angle));
        if(tmp<neardis)
        {
            n=s[i];
            neardis=tmp;
        }
    }
    return n;
}

FILE* fino;
FILE* findp;
FILE* fincnf;
FILE* fout;
double crange;
double rrange;
double csize;
point center;
bool ca[1000][1000];
int numobs,numdps,numk,numcells;
double osa,osd;
vector<point> dpsdata;
vector<obs> obsdata;
vector<triangle> tempd;
vector< vector<triangle> > vistri;
vector< vector<triangle> > res;
vector< vector<triangle> > clists;
vector< vector<int> > al;
vector<triangle> oldlist;
triangle a,b;

void ps(seg s)
{
    printf("%d %lf %lf %lf %lf\n",s.id,s.p1.x,s.p1.y,s.p2.x,s.p2.y);
}

//trivial rejection korte hobe

vector<seg> gentri(vector<obs> o)
{
    vector<seg> ret,sl,sids;
    vector<int> tmp;
    int i,j,sid=0;
    //for(i=0;i<o.size();i++)cout<<o[i].xl<<endl;
    for(i=0;i<o.size();i++)
    {
        if(fabs(o[i].xl)<eps)
        {
            sl.push_back(seg(point(o[i].xh,o[i].yl),point(0,o[i].yl),true,true,sid));
            sids.push_back(seg(point(o[i].xh,o[i].yl),point(0,o[i].yl),true,true,sid));
            sl.push_back(seg(point(0,o[i].yl),point(o[i].xh,o[i].yl),true,false,sid));
            sid++;
        }
        else if(fabs(o[i].yl)<eps)
        {

            sl.push_back(seg(point(o[i].xl,0),point(o[i].xl,o[i].yh),false,true,sid));
            sids.push_back(seg(point(o[i].xl,0),point(o[i].xl,o[i].yh),false,true,sid));
            sl.push_back(seg(point(o[i].xl,o[i].yh),point(o[i].xl,0),false,false,sid));
            sid++;
        }
        else
        {
            sl.push_back(seg(point(o[i].xh,o[i].yl),point(o[i].xl,o[i].yl),true,true,sid));
            sids.push_back(seg(point(o[i].xh,o[i].yl),point(o[i].xl,o[i].yl),true,true,sid));
            sl.push_back(seg(point(o[i].xl,o[i].yl),point(o[i].xh,o[i].yl),true,false,sid));
            sid++;

            sl.push_back(seg(point(o[i].xl,o[i].yl),point(o[i].xl,o[i].yh),false,true,sid));
            sids.push_back(seg(point(o[i].xl,o[i].yl),point(o[i].xl,o[i].yh),false,true,sid));
            sl.push_back(seg(point(o[i].xl,o[i].yh),point(o[i].xl,o[i].yl),false,false,sid));
            sid++;
        }
    }
    sl.push_back(seg(point(crange+eps,0),point(crange+eps,crange+eps),false,true,sid));
    sids.push_back(seg(point(crange+eps,0),point(crange+eps,crange+eps),false,true,sid));
    sl.push_back(seg(point(crange+eps,crange+eps),point(crange+eps,0),false,false,sid));
    sid++;

    sl.push_back(seg(point(crange+eps,crange+eps),point(0,crange+eps),true,true,sid));
    sids.push_back(seg(point(crange+eps,crange+eps),point(0,crange+eps),true,true,sid));
    sl.push_back(seg(point(0,crange+eps),point(crange+eps,crange+eps),true,false,sid));
    sid++;


    sort(sl.begin(),sl.end(),cmpseg);
    //for(i=0;i<sl.size();i++)ps(sl[i]);
    //for(i=0;i<sids.size();i++)ps(sids[i]);

    double ang=angle1(sl[0]);
    tmp.clear();
    tmp.push_back(sl[0].id);

    for(i=1;i<sl.size();i++)
    {
        //cout<<sl.size()<<endl;
        double nang=angle1(sl[i]);
        //cout<<nang<<endl;
        if(fabs(nang-ang)>eps)
        {
            vector<seg> param;
            param.clear();
            for(j=0;j<tmp.size();j++)param.push_back(sids[tmp[j]]);
            seg det;
            if(param.size()==0)det=sl[i];
            else det=nearinter(param,nang);
            ret.push_back(seg(intersect(det,ang),intersect(det,nang)));
        }
        ang=nang;
        if(sl[i].dir==true)
        {
            tmp.push_back(sl[i].id);
        }
        else
        {
            int k;
            for(k=0;k<tmp.size();k++)if(tmp[k]==sl[i].id)break;
            tmp.erase(tmp.begin()+k);
        }
    }
    return ret;
}

vector<triangle> genvisreg(point c,vector<obs> ol)
{
    vector<triangle> ret;
    vector<obs> ol1,ol2,ol3,ol4,tmpo;
    vector<seg> s1,s2,s3,s4,tmps;
    int i;
    for(i=0;i<ol.size();i++)
    {
        ol[i].xl-=c.x;
        if(ol[i].xl<-crange)ol[i].xl=-crange;
        ol[i].xh-=c.x;
        if(ol[i].xh>crange)ol[i].xh=crange;
        ol[i].yl-=c.y;
        if(ol[i].yl<-crange)ol[i].yl=-crange;
        ol[i].yh-=c.y;
        if(ol[i].yh>crange)ol[i].yh=crange;
    }

    for(i=0;i<ol.size();i++)
    {
        if     (ol[i].xl>=0&&ol[i].yl>=0)ol1.push_back(ol[i]);
        else if(ol[i].xh<=0&&ol[i].yl>=0)ol2.push_back(ol[i]);
        else if(ol[i].xh<=0&&ol[i].yh<=0)ol3.push_back(ol[i]);
        else if(ol[i].xl>=0&&ol[i].yh<=0)ol4.push_back(ol[i]);
        else if(ol[i].xl<=0&&ol[i].xh>=0)
        {
            if(ol[i].yl>=0)
            {
                ol1.push_back(obs(0,ol[i].xh,ol[i].yl,ol[i].yh));
                ol2.push_back(obs(ol[i].xl,0,ol[i].yl,ol[i].yh));
            }
            else
            {
                ol4.push_back(obs(0,ol[i].xh,ol[i].yl,ol[i].yh));
                ol3.push_back(obs(ol[i].xl,0,ol[i].yl,ol[i].yh));
            }
        }
        else if(ol[i].yl<=0&&ol[i].yh>=0)
        {
            if(ol[i].xl>=0)
            {
                ol1.push_back(obs(ol[i].xl,ol[i].xh,0,ol[i].yh));
                ol4.push_back(obs(ol[i].xl,ol[i].xh,ol[i].yl,0));
            }
            else
            {
                ol2.push_back(obs(ol[i].xl,ol[i].xh,0,ol[i].yh));
                ol3.push_back(obs(ol[i].xl,ol[i].xh,ol[i].yl,0));
            }
        }

    }

    s1=gentri(ol1);

    //cout<<s1.size()<<endl;

    tmpo.clear();
    tmps.clear();
    for(i=0;i<ol2.size();i++)tmpo.push_back(obsrot(ol2[i],-pi/2.0));
    tmps=gentri(tmpo);
    s2.clear();
    for(i=0;i<tmps.size();i++)s2.push_back(segrot(tmps[i],pi/2.0));

    tmpo.clear();
    tmps.clear();
    for(i=0;i<ol3.size();i++)tmpo.push_back(obsrot(ol3[i],-pi));
    tmps=gentri(tmpo);
    s3.clear();
    for(i=0;i<tmps.size();i++)s3.push_back(segrot(tmps[i],pi));

    tmpo.clear();
    tmps.clear();
    for(i=0;i<ol4.size();i++)tmpo.push_back(obsrot(ol4[i],pi/2.0));
    tmps=gentri(tmpo);
    s4.clear();
    for(i=0;i<tmps.size();i++)s4.push_back(segrot(tmps[i],-pi/2.0));

    ret.clear();
    for(i=0;i<s1.size();i++)ret.push_back(triangle(point(c.x,c.y),point(c.x+s1[i].p1.x,c.y+s1[i].p1.y),point(c.x+s1[i].p2.x,c.y+s1[i].p2.y)));
    for(i=0;i<s2.size();i++)ret.push_back(triangle(point(c.x,c.y),point(c.x+s2[i].p1.x,c.y+s2[i].p1.y),point(c.x+s2[i].p2.x,c.y+s2[i].p2.y)));
    for(i=0;i<s3.size();i++)ret.push_back(triangle(point(c.x,c.y),point(c.x+s3[i].p1.x,c.y+s3[i].p1.y),point(c.x+s3[i].p2.x,c.y+s3[i].p2.y)));
    for(i=0;i<s4.size();i++)ret.push_back(triangle(point(c.x,c.y),point(c.x+s4[i].p1.x,c.y+s4[i].p1.y),point(c.x+s4[i].p2.x,c.y+s4[i].p2.y)));


    return ret;
}


void readfilet()
{
    fino=fopen("datat.txt","r");

	if(fino==NULL)
		printf("Error opening file");

	fscanf(fino," %lf %lf %lf %lf %lf %lf",&a.a.x,&a.a.y,&a.b.x,&a.b.y,&a.c.x,&a.c.y);
	fscanf(fino," %lf %lf %lf %lf %lf %lf",&b.a.x,&b.a.y,&b.b.x,&b.b.y,&b.c.x,&b.c.y);

    fclose(fino);
}

void readfile()
{

    int i;
    point p;
    double xl,xh,yl,yh;



    fincnf=fopen("config.txt","r");

	if(fincnf==NULL)
		printf("Error opening file");

	fscanf(fincnf," %lf",&rrange);
	fscanf(fincnf," %lf",&crange);
    fscanf(fincnf," %d %d",&numdps,&numk);
	fscanf(fincnf," %d",&numobs);
	fscanf(fincnf," %lf %lf",&osa,&osd);
	fscanf(fincnf," %d",&numcells);

	csize=rrange/(double)numcells;

    fclose(fincnf);



    fino=fopen("dataobs.txt","r");

	if(fino==NULL)
		printf("Error opening file");

    obsdata.clear();
	for(i=0;i<numobs;i++)
    {
        fscanf(fino," %lf %lf %lf %lf",&xl,&xh,&yl,&yh);
        obsdata.push_back(obs(xl,xh,yl,yh));
    }
    fclose(fino);



    findp=fopen("datadps.txt","r");

	if(findp==NULL)
		printf("Error opening file");

	dpsdata.clear();
	for(i=0;i<numdps;i++)
    {
        fscanf(findp," %lf %lf",&p.x,&p.y);
        dpsdata.push_back(p);
    }

    fclose(findp);


}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
			break;

		default:
			break;
	}
}

void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			break;
		case GLUT_KEY_UP:		// up arrow key
			break;

		case GLUT_KEY_RIGHT:
			break;
		case GLUT_KEY_LEFT:
			break;

		case GLUT_KEY_PAGE_UP:
			break;
		case GLUT_KEY_PAGE_DOWN:
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:

			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(rrange/2.0,rrange/2.0,rrange/2.0,	rrange/2.0,rrange/2.0,0,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	int i,j;

    glColor3f(1,0,0);
/*
    for(i=0;i<vistri.size();i++)
    {
        for(j=0;j<vistri[i].size();j++)
        {
            glBegin(GL_LINES);
            {
                glVertex3f(vistri[i][j].a.x,vistri[i][j].a.y,0);
                glVertex3f(vistri[i][j].b.x,vistri[i][j].b.y,0);
                glVertex3f(vistri[i][j].b.x,vistri[i][j].b.y,0);
                glVertex3f(vistri[i][j].c.x,vistri[i][j].c.y,0);
                glVertex3f(vistri[i][j].c.x,vistri[i][j].c.y,0);
                glVertex3f(vistri[i][j].a.x,vistri[i][j].a.y,0);

            }
            glEnd();

        }

    }
*/
    glColor3f(0,1,0);
    for(i=0;i<tempd.size();i++)
    {
        glBegin(GL_LINES);
        {
            glVertex3f(tempd[i].a.x,tempd[i].a.y,0);
            glVertex3f(tempd[i].b.x,tempd[i].b.y,0);
            glVertex3f(tempd[i].b.x,tempd[i].b.y,0);
            glVertex3f(tempd[i].c.x,tempd[i].c.y,0);
            glVertex3f(tempd[i].c.x,tempd[i].c.y,0);
            glVertex3f(tempd[i].a.x,tempd[i].a.y,0);

        }
        glEnd();

    }

    glColor3f(0,0,1);
    for(i=0;i<dpsdata.size();i++)
    {
        glBegin(GL_QUADS);
        {
			glVertex3f(dpsdata[i].x+5,dpsdata[i].y+5,0);
			glVertex3f(dpsdata[i].x-5,dpsdata[i].y+5,0);
			glVertex3f(dpsdata[i].x-5,dpsdata[i].y-5,0);
			glVertex3f(dpsdata[i].x+5,dpsdata[i].y-5,0);
        }
        glEnd();
    }

    glColor3f(1,0,0);
    for(i=0;i<obsdata.size();i++)
    {
        glBegin(GL_QUADS);
        {
			glVertex3f(obsdata[i].xl,obsdata[i].yl,0);
			glVertex3f(obsdata[i].xh,obsdata[i].yl,0);
			glVertex3f(obsdata[i].xh,obsdata[i].yh,0);
			glVertex3f(obsdata[i].xl,obsdata[i].yh,0);
        }
        glEnd();
    }

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}


bool doesint1(obs a,obs b)
{
    if(((b.xl+eps<a.xl&&a.xl+eps<b.xh)||(b.xl+eps<a.xh&&a.xh+eps<b.xh)||(a.xl+eps<b.xl&&b.xl+eps<a.xh)||(a.xl+eps<b.xh&&b.xh+eps<a.xh))&&((b.yl+eps<a.yl&&a.yl+eps<b.yh)||(b.yl+eps<a.yh&&a.yh+eps<b.yh)||(a.yl+eps<b.yl&&b.yl+eps<a.yh)||(a.yl+eps<b.yh&&b.yh+eps<a.yh)))return true;

    return false;
}

vector<obs> rangequery(point p)
{
    int j;
    vector<obs> ret;
    obs to;
    to.xl=p.x-crange;
    to.xh=p.x+crange;
    to.yl=p.y-crange;
    to.yh=p.y+crange;
    for(j=0;j<obsdata.size();j++)
    {
        if(doesint1(obsdata[j],to)==true)ret.push_back(obsdata[j]);
    }
    return ret;
}

pair<bool,double> scinter(double sc,point a,point b)
{
    if(a.y<b.y)
    {
        point tmp=a;
        a=b;
        b=tmp;
    }
    if(sc>a.y||sc<b.y)return pair<bool,double>(false,0.0);
    point ret=intp(a,b,point(0,sc),point(rrange,sc));
    return pair<bool,double>(true,ret.x);
}

int aaa[10000];

double countarea(unsigned int bm)
{
    int i,j;
    double ret=0.0;

    for(i=0;i<oldlist.size();i++)
    {
        triangle ttt=oldlist[i];
        if(((~bm)&ttt.bitmap)==0)
        {
            aaa[i]++;
            ret+=ttt.area;
        }
    }
    return ret;
}

vector<double> subset(vector<int> s)
{
    int i,j;
    unsigned int mask;
    int sz=s.size();
    vector< pair<unsigned int,double> > mc;
    mc.clear();
    for(i=1;i<(1<<sz);i++)
    {
        mask=0;
        for(j=0;j<sz;j++)
        {
            if(((1<<j)&i)!=0)mask=(mask|(1<<s[j]));
        }
        mc.push_back(pair<unsigned int,double>(mask,countarea(mask)));
    }
    vector<double> ret;
    ret.clear();
    for(i=0;i<=sz;i++)ret.push_back(-1);
    for(i=0;i<mc.size();i++)
    {
        int none=numone(mc[i].first);
        if(ret[none]<mc[i].second)ret[none]=mc[i].second;
    }
    ret[0]=0;

    vector<double> rett;
    rett.clear();
    for(i=0;i<ret.size()-1;i++)rett.push_back(ret[i+1]-ret[i]);
    return rett;
}


bool marked[32];

vector<int> dfs(int n)
{
    vector<int> ret;
    ret.clear();
    ret.push_back(n);
    int i,j;
    marked[n]=true;
    for(i=0;i<al[n].size();i++)
    {
        vector<int> tempp;
        tempp.clear();
        if(marked[al[n][i]]==false)tempp=dfs(al[n][i]);
        for(j=0;j<tempp.size();j++)ret.push_back(tempp[j]);
    }
    return ret;
}

vector< vector<int> > cc()
{
    int i;
    vector< vector<int> > ret;
    for(i=0;i<numdps;i++)marked[i]=false;
    for(i=0;i<al.size();i++)if(marked[i]==false)ret.push_back(dfs(i));
    return ret;
}

unsigned int revbit(unsigned int n)
{
    unsigned int ret=0;
    int i;
    for(i=0;i<numdps;i++)
    {
        if(n%2==1)ret=ret+(1<<(numdps-1-i));
        n/=2;
    }
    return ret;
}


void init(){
	//codes for initialization
    int i,j,k,l;
    readfile();
    fout=fopen("outc.txt","w+");
    if(fout==NULL)
		printf("Error opening file");


    const clock_t begin_time = clock();

    for(i=0;i<numdps;i++)vistri.push_back(genvisreg(dpsdata[i],rangequery(dpsdata[i])));


    oldlist.clear();
    vector<triangle> common,onotn,nnoto;
    oldlist=vistri[0];
    for(i=0;i<oldlist.size();i++)
    {
        oldlist[i].bitmap=1;
        oldlist[i].cnt=1;
    }
    for(i=1;i<numdps;i++)
    {
        common.clear();
        for(j=0;j<vistri[i].size();j++)
        {
            for(k=0;k<oldlist.size();k++)
            {
                triangle tn,to;
                tn=vistri[i][j];
                to=oldlist[k];
                res=inttri(to,tn);
                triangle ret;
                for(l=0;l<res[0].size();l++)
                {
                    ret=res[0][l];
                    ret.bitmap=2*to.bitmap+1;
                    ret.cnt=to.cnt+1;
                    common.push_back(ret);
                }
            }
        }



        onotn=oldlist;
        for(j=0;j<vistri[i].size();j++)
        {
            triangle tn,to;
            tn=vistri[i][j];
            vector<triangle> tmp=onotn;
            onotn.clear();
            for(k=0;k<tmp.size();k++)
            {
                to=tmp[k];
                res=inttri(to,tn);
                triangle ret;
                for(l=0;l<res[1].size();l++)
                {
                    ret=res[1][l];
                    ret.bitmap=to.bitmap;
                    ret.cnt=to.cnt;
                    onotn.push_back(ret);
                }
            }
        }
        for(l=0;l<onotn.size();l++)onotn[l].bitmap=2*onotn[l].bitmap;




        nnoto=vistri[i];
        for(j=0;j<oldlist.size();j++)
        {
            vector<triangle> tmp=nnoto;
            nnoto.clear();
            for(k=0;k<tmp.size();k++)
            {
                triangle tn,to;
                tn=oldlist[j];
                to=tmp[k];
                res=inttri(tn,to);
                triangle ret;
                for(l=0;l<res[2].size();l++)
                {
                    ret=res[2][l];
                    nnoto.push_back(ret);
                }
            }

        }
        for(l=0;l<nnoto.size();l++)nnoto[l].bitmap=nnoto[l].cnt=1;

        oldlist.clear();
        for(j=0;j<common.size();j++)oldlist.push_back(common[j]);
        for(j=0;j<onotn.size();j++)oldlist.push_back(onotn[j]);
        for(j=0;j<nnoto.size();j++)oldlist.push_back(nnoto[j]);

    }
    tempd.clear();
    tempd=oldlist;
    for(i=0;i<oldlist.size();i++)
    {
        oldlist[i].bitmap=revbit(oldlist[i].bitmap);
        oldlist[i].area=0.5*(fabs(area3p(oldlist[i].a,oldlist[i].b,oldlist[i].c)));
    }

    map<unsigned int,double> mp;
    map<unsigned int,double>::iterator it;
    for(i=0;i<oldlist.size();i++)
    {
        it=mp.find(oldlist[i].bitmap);
        if(it==mp.end())
        {
            mp.insert(pair<unsigned int,double>(oldlist[i].bitmap,0.5*fabs(areat(oldlist[i]))));
        }
        else mp[oldlist[i].bitmap]+=0.5*fabs(areat(oldlist[i]));
    }

    vector<record> recs;
    for(it=mp.begin();it!=mp.end();it++)
    {
        unsigned int bm=it->first;
        double areaa=it->second;
        recs.push_back(record(bm,areaa,numone(bm)));
        //printf("%u %lf\n",bm,areaa);
    }

    for(i=0;i<numdps;i++)
    {
        vector<triangle> tempp;
        tempp.clear();
        clists.push_back(tempp);
    }



    bool ag[32][32];
    for(i=0;i<32;i++)for(j=0;j<32;j++)ag[i][j]=false;

    for(i=0;i<numdps;i++)
    {
        for(j=i+1;j<numdps;j++)
        {
            if(i!=j)
            {
                for(k=0;k<vistri[i].size();k++)
                {
                    for(l=0;l<vistri[j].size();l++)
                    {
                        triangle ti=vistri[i][k];
                        triangle tj=vistri[j][l];
                        res=inttri(ti,tj);
                        if(res[0].size()>0)
                        {
                            ag[i][j]=true;
                            ag[j][i]=true;
                            break;
                        }
                    }
                }
            }
        }
    }


    al.clear();
    for(i=0;i<numdps;i++)
    {
        vector<int> tmpp;
        tmpp.clear();
        for(j=0;j<numdps;j++)if(ag[i][j]==true)tmpp.push_back(j);
        al.push_back(tmpp);
    }




    vector< vector<int> > comps;
    comps.clear();
    comps=cc();


//    printf("\n");
//    for(i=0;i<comps.size();i++)
//    {
//        for(j=0;j<comps[i].size();j++)printf("%d ",comps[i][j]);
//        printf("\n");
//    }



    vector<double> aa,bb;
    aa.clear();
    int numc=comps.size();
    for(i=0;i<numc;i++)
    {
        bb=subset(comps[i]);
        for(j=0;j<bb.size();j++)aa.push_back(bb[j]);
    }
    sort(aa.begin(),aa.end());
    reverse(aa.begin(),aa.end());

    const clock_t end_time = clock();

    fprintf(fout,"%d\n",numdps);
    for(i=1;i<aa.size();i++)aa[i]=aa[i]+aa[i-1];
    for(i=0;i<aa.size();i++)fprintf(fout,"%.5lf\n",aa[i]);


    cout<<"Time Difference: "<<end_time-begin_time<<" miliseconds"<<endl;


    fclose(fout);

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(90,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(676, 676);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
