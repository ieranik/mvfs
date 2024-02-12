#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <map>
#include <conio.h>

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

vector<obs> obsdata;
vector<point> dpsdata;
vector<obs> tmp;
int numobs,rrange,crange,osa,osd,numdps,numcells,numk;
FILE* fin;
FILE *fout;
FILE *fout1;

void po(obs o)
{
    cout<<o.xl<<" "<<o.xh<<" "<<o.yl<<" "<<o.yh<<endl;
}

bool doesint1(obs a,obs b)
{
    if(((b.xl+eps<a.xl&&a.xl+eps<b.xh)||(b.xl+eps<a.xh&&a.xh+eps<b.xh)||(a.xl+eps<b.xl&&b.xl+eps<a.xh)||(a.xl+eps<b.xh&&b.xh+eps<a.xh))&&((b.yl+eps<a.yl&&a.yl+eps<b.yh)||(b.yl+eps<a.yh&&a.yh+eps<b.yh)||(a.yl+eps<b.yl&&b.yl+eps<a.yh)||(a.yl+eps<b.yh&&b.yh+eps<a.yh)))return true;

    return false;
}

bool doesintn(vector<obs> ol,obs o)
{
    int i;
    for(i=0;i<ol.size();i++)if(doesint1(ol[i],o)==true)return true;

    return false;
}

bool isinside1(obs o,point p)
{
    if(o.xl<p.x+eps&&p.x<o.xh+eps&&o.yl<p.y+eps&&p.y<o.yh+eps)return true;

    return false;
}

bool isinsiden(vector<obs> ol,point p)
{
    int i;
    for(i=0;i<ol.size();i++)if(isinside1(ol[i],p)==true)return true;

    return false;
}

void genobs()
{
    obsdata.clear();
    fout=fopen("dataobs.txt","w+");
    if(fout==NULL)
		printf("Error opening file");

    obs o;
    int i,ri,rii;
    double rd;
    srand(time(NULL));
    int cnt=numobs;
    while(cnt!=0)
    {
        ri=rand()%(rrange-osa-osd);
        rii=rand()%1000;
        o.xl=ri+(double)rii/1000.0;
        ri=rand()%(2*osd)+(osa-osd);
        rii=rand()%1000;
        o.xh=o.xl+ri+(double)rii/1000.0;

        ri=rand()%(rrange-osa-osd);
        rii=rand()%1000;
        o.yl=ri+(double)rii/1000.0;
        ri=rand()%(2*osd)+(osa-osd);
        rii=rand()%1000;
        o.yh=o.yl+ri+(double)rii/1000.0;


        if(doesintn(obsdata,o)==false)
        {
            fprintf(fout,"%.3lf\t%.3lf\t%.3lf\t%.3lf\n",o.xl,o.xh,o.yl,o.yh);
            obsdata.push_back(o);
            cnt--;
            //cout<<cnt<<endl;
        }
    }
    fclose(fout);
}

void gendps()
{
    dpsdata.clear();
    fout1=fopen("datadps.txt","w+");
    if(fout1==NULL)
		printf("Error opening file");

    int i,ri,rii;
    point p;
    int cnt=numdps;
    while(cnt!=0)
    {
        ri=rand()%(rrange-2*crange);
        rii=rand()%1000;
        p.x=crange+ri+(double)rii/1000.0;

        ri=rand()%(rrange-2*crange);
        rii=rand()%1000;
        p.y=crange+ri+(double)rii/1000.0;

        if(isinsiden(obsdata,p)==false)
        {
            fprintf(fout1,"%.3lf\t%.3lf\n",p.x,p.y);
            dpsdata.push_back(p);
            cnt--;
        }

    }
    fclose(fout1);
}

void readfile()
{

    fin=fopen("config.txt","r");

	if(fin==NULL)
		printf("Error opening file");

	fscanf(fin," %d",&rrange);
	fscanf(fin," %d",&crange);

    fscanf(fin," %d %d",&numdps,&numk);

	fscanf(fin," %d",&numobs);
    fscanf(fin," %d %d",&osa,&osd);

    fscanf(fin," %d",&numcells);

    fclose(fin);

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
	gluLookAt(rrange/2,rrange/2,rrange/2,	rrange/2,rrange/2,0,	0,1,0);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	int i,j;

    glColor3f(0,1,0);

    for(i=0;i<tmp.size();i++)
    {
        glBegin(GL_QUADS);
        {
			glVertex3f(tmp[i].xl,tmp[i].yl,0);
			glVertex3f(tmp[i].xh,tmp[i].yl,0);
			glVertex3f(tmp[i].xh,tmp[i].yh,0);
			glVertex3f(tmp[i].xl,tmp[i].yh,0);
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

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
    int i,j,k,l;

    readfile();
    genobs();
    gendps();



    for(i=0;i<dpsdata.size();i++)
    {
        obs to;
        to.xl=dpsdata[i].x-crange;
        to.xh=dpsdata[i].x+crange;
        to.yl=dpsdata[i].y-crange;
        to.yh=dpsdata[i].y+crange;
        for(j=0;j<obsdata.size();j++)
        {
            if(doesint1(obsdata[j],to)==true)tmp.push_back(obsdata[j]);
        }
    }
    //cout<<tmp.size();

    //getch();




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
	gluPerspective(90,	1,	1,	10500.0);
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
