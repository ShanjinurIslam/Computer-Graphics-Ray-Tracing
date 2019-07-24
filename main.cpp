//
//  main.cpp
//  Computer-Graphics-Ray-Tracing
//
//  Created by Shanjinur Islam on 7/17/19.
//  Copyright © 2019 Shanjinur Islam. All rights reserved.
//
#define GL_SILENCE_DEPRECATION
#include <iostream>
#include <windows.h>
#include <glut.h>
#include <cmath>
#include <fstream>
#include <vector>
#include "bitmap_image.hpp"

#define GL_SILENCE_DEPRECATION
#define pi acos(-1.0)

using namespace std ;

double width, height;
double recursionLevel;
double cameraHeight;
double cameraAngle;
double angle;
double increment;
int drawgrid;
int drawaxes;

class Point{
public:
    double x,y,z ;
    Point(double x=0,double y=0,double z=0){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Point operator+(Point v)
    {
        Point temp;
        temp.x = this->x + v.x;
        temp.y = this->y + v.y;
        temp.z = this->z + v.z;
        return temp;
    }

    Point &operator=(Point v)
    {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;

        return *this;
    }

}pos;

class Vector
{
public:
    double ax, ay, az;

    Vector(double ax=0, double ay=0, double az=0)
    {
        this->ax = ax;
        this->ay = ay;
        this->az = az;
    }

    void normalize()
    {
        double r = sqrt((ax * ax) + (ay * ay) + (az * az));
        ax = ax / r;
        ay = ay / r;
        az = az / r;
    }

    Vector operator+(Vector v)
    {
        Vector temp;
        temp.ax = this->ax + v.ax;
        temp.ay = this->ay + v.ay;
        temp.az = this->az + v.az;

        return temp;
    }

    Vector operator-(Vector v)
    {
        Vector temp;
        temp.ax = this->ax - v.ax;
        temp.ay = this->ay - v.ay;
        temp.az = this->az - v.az;

        return temp;
    }

    Vector &operator=(Vector v)
    {
        this->ax = v.ax;
        this->ay = v.ay;
        this->az = v.az;

        return *this;
    }

    Vector operator*(double d)
    {
        Vector temp;
        temp.ax = ax * d;
        temp.ay = ay * d;
        temp.az = az * d;
        return temp;
    }
    void print()
    {
        cout << ax << " " << ay << " " << az << endl;
    }

    double dot(Vector x)
    {
        return (this->ax * x.ax) + (this->ay * x.ay) + (this->az * x.az);
    }

    Vector cross(Vector x)
    {
        Vector temp;
        temp.ax = (ay * x.az - az * x.ay);
        temp.ay = (az * x.ax - ax * x.az);
        temp.az = (ax * x.ay - ay * x.ax);
        return temp;
    }

    Vector generateVector(Point a,Point b){ //to generate vector AB
        Vector temp ;
        temp.ax = b.x - a.x ;
        temp.ay = b.y - a.y ;
        temp.az = b.z - a.z ;

        return temp ;
    }
}u,l,r;

class Ray{
    Point point ;
    Vector vector ;
    Ray(){
        this->point = Point() ;
        this->vector = Vector() ;
    }
    Ray(Point p,Vector v){
        this->point = p ;
        this->vector = v ;
    }

    Ray(Point a,Point b){
        vector = vector.generateVector(a,b) ;
    }
};

class Color{
public:
    double red,green,blue ;
    Color(double red=0,double green=0,double blue=0){
        this->red   = red ;
        this->green = green ;
        this->blue  = blue ;
    }
};

class Triangle{
public:
    Point a,b,c ;
    Color color;

    Triangle(){
    }

    Triangle(Point a,Point b,Point c,Color color){
        this->a = a ;
        this->b = b ;
        this->c = c ;
        this->color = color ;
    }

    Triangle &operator=(Triangle v)
    {
        this->a = v.a ;
        this->b = v.b ;
        this->c = v.c ;
        this->color = v.color ;

        return *this;
    }

    void draw(){
        glColor3d(color.red,color.green,color.blue) ;
        glBegin(GL_TRIANGLES);
        glVertex3d(a.x,a.y,a.z);
        glVertex3d(b.x,b.y,b.z);
        glVertex3d(c.x,c.y,c.z);
        glEnd();
    }
};

class Square{
public:
    Point *corner ;
    Color color;
    Square(){
        corner = new Point[4] ;
        for(int i=0;i<4;i++){
            corner[i] = Point() ;
        }
        color = Color() ;
    }

    Square(Point start,double side_length,Color color){
        corner = new Point[4] ;
        for(int i=0;i<4;i++){
            corner[i] = Point() ;
        }
        this->color = color ;
        corner[0] = start ;
        corner[1] = corner[0] + Point(side_length,0,0);
        corner[2] = corner[1] + Point(0,side_length,0);
        corner[3] = corner[0] + Point(0,side_length,0);
    }

    Square &operator=(Square s)
    {
        for(int i=0;i<4;i++){
            corner[i] = s.corner[i] ;
        }
        color = s.color ;
        return *this ;
    }

    void draw() {
        glColor3f(color.red,color.green,color.blue);
        glBegin(GL_QUADS);
        glVertex3f(corner[0].x,corner[0].y,corner[0].z);
        glVertex3f(corner[1].x,corner[1].y,corner[1].z);
        glVertex3f(corner[2].x,corner[2].y,corner[2].z);
        glVertex3f(corner[3].x,corner[3].y,corner[3].z);
        glEnd();
    }
    Point* getCornerPoints(){
        return corner ;
    }
};

class Pyramid{
public:
    Point low,top;
    double height,width ;
    Color color ;
    Triangle sides[4];
    Square bottom ;

    Pyramid(Point low,double width,double height,Color color){
        this->low = low ;
        this->height = height ;
        this->width = width ;

        for(int i=0;i<4;i++){
            sides[i] = Triangle() ;
        }

        top = low + Point(width/2,width/2,height);
        bottom = Square(low,width,color) ;

        Point *corners = new Point[4] ;
        corners = bottom.getCornerPoints() ;

        for(int i=0;i<4;i++){
            sides[i] = Triangle(corners[i],corners[(i+1)%4],top,color) ;
        }
    }

    void draw(){
        bottom.draw() ;
        for(int i=0;i<4;i++){
            sides[i].draw();
        }
    }
};

class Sphere{
public:
    Point center ;
    double radius;
    Color color ;

    Sphere(Point center,double radius,Color color){
        this->center = center ;
        this->radius = radius ;
        this->color = color ;
    }

    void drawSphere(double radius, int slices, int stacks)
    {
        Point p[stacks + 1][slices + 1];
        double h, r, angle1, angle2;
        for (int i = 0; i <= stacks; i++)
        {
            angle1 = ((double)i / (double)stacks) * (acos(-1.0)*2);
            h = radius * sin(angle1);
            r = radius * cos(angle1);
            for (int j = 0; j <= slices; j++)
            {
                angle2 = ((double)j / (double)slices) * acos(-1.0) * 2;
                p[i][j].x = r * cos(angle2);
                p[i][j].y = r * sin(angle2);
                p[i][j].z = h  ;
            }
        }
        for (int i = 0; i < stacks; i++)
        {
            for (int j = 0; j < slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    glVertex3f(p[i][j].x, p[i][j].y, p[i][j].z);
                    glVertex3f(p[i][j + 1].x, p[i][j + 1].y, p[i][j + 1].z);
                    glVertex3f(p[i + 1][j + 1].x, p[i + 1][j + 1].y, p[i + 1][j + 1].z);
                    glVertex3f(p[i + 1][j].x, p[i + 1][j].y, p[i + 1][j].z);
                }
                glEnd();
            }
        }
        glPopMatrix();
    }

    void draw(){
        glPushMatrix() ;
        glTranslated(center.x,center.y,center.z);
        glColor3f(color.red,color.green,color.blue);
        drawSphere(radius,25,25);
        glPopMatrix() ;
    }
};

class ChessBoard{
private:
    int color_map[101][101];
public:
    ChessBoard(){
        int flag = 1 ;
        for(int i=0;i<=100;i++){
            for(int j=0;j<=100;j++){
                color_map[i][j] = flag%2 ;
                flag++;
            }
        }
    }
    void draw(){
        int flag=1;
        for(double i=-30*50;i<=30*50;i+=30){
            for(double j=-30*50;j<=30*50;j+=30){
                if(color_map[(int)(i/30.0)+50][(int)(j/30.0)+50]){
                    Square(Point(i,j,0),30,Color(0,0,0)).draw();
                }
                else{
                    Square(Point(i,j,0),30,Color(255,255,255)).draw();
                }
                flag++;
            }
        }
    }
};

vector<Pyramid> pyramids ;
vector<Sphere> spheres ;
ChessBoard chessboard ;


//test purpose
void drawAxes()
{
    glPushMatrix();
    glBegin(GL_LINES);
    {
        glColor3f(1.0, 0, 0);
        glVertex3f(100, 0, 0);
        glVertex3f(-100, 0, 0);

        glColor3f(0, 1.0, 0);
        glVertex3f(0, -100, 0);
        glVertex3f(0, 100, 0);

        glColor3f(0, 0, 1.0);
        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, -100);
    }
    glEnd();
    glPopMatrix();
}

void init()
{
    pos = Point(50, 50, 50);
    u = Vector(0, 0, 1); //z axis is up vector
    l = Vector(-1 / (sqrt(2)), -1 / sqrt(2), 0);
    r = Vector(-1 / (sqrt(2)), 1 / sqrt(2), 0);
    angle = acos(-1.0) / 30; //3 degree angle change
    increment = 5 ;
    chessboard = ChessBoard();
    glClearColor(0, 0, 0, 0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(80, 1, 1, 1000.0);
}


void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(pos.x, pos.y, pos.z, pos.x + l.ax, pos.y + l.ay, pos.z + l.az, u.ax, u.ay, u.az);

    glMatrixMode(GL_MODELVIEW);

    glPushMatrix();
    chessboard.draw() ;

    for(int i=0;i<pyramids.size();i++){
        pyramids[i].draw() ;
    }

    for(int i=0;i<spheres.size();i++){
        spheres[i].draw() ;
    }

    glPopMatrix();

    glutSwapBuffers();
}
void animate()
{
    glutPostRedisplay();
}
void keyboardListener(unsigned char key, int x, int y)
{
    /*
     U
     |
     |
     |
     L ______ R

     1. Rotating will not change the UP vector but will change the LOOK Vector and RIGHT Vector.
     2. Looking up/down not change the RIGHT vector but will change the LOOK Vector and RIGHT Vector.
     3. Tilt not change the LOOK vector but will change the UP Vector and RIGHT Vector.
     */

    switch (key)
    {
        case '1':
            r.ax = r.ax * cos(-1.0 * angle) + l.ax * sin(-1.0 * angle);
            r.ay = r.ay * cos(-1.0 * angle) + l.ay * sin(-1.0 * angle);
            r.az = r.az * cos(-1.0 * angle) + l.az * sin(-1.0 * angle);
            l = u.cross(r) ;
            break;
        case '2':
            r.ax = r.ax * cos(angle) + l.ax * sin(angle);
            r.ay = r.ay * cos(angle) + l.ay * sin(angle);
            r.az = r.az * cos(angle) + l.az * sin(angle);
            l = u.cross(r);
            break;
        case '3':
            l.ax = l.ax * cos(angle) + u.ax * sin(angle);
            l.ay = l.ay * cos(angle) + u.ay * sin(angle);
            l.az = l.az * cos(angle) + u.az * sin(angle);
            u = r.cross(l) ;
            break;
        case '4':
            l.ax = l.ax * cos(-1.0 * angle) + u.ax * sin(-1.0 * angle);
            l.ay = l.ay * cos(-1.0 * angle) + u.ay * sin(-1.0 * angle);
            l.az = l.az * cos(-1.0 * angle) + u.az * sin(-1.0 * angle);
            u = r.cross(l) ;
            break;
        case '5':
            u.ax = u.ax * cos(angle) + r.ax * sin(angle);
            u.ay = u.ay * cos(angle) + r.ay * sin(angle);
            u.az = u.az * cos(angle) + r.az * sin(angle);
            r = l.cross(u);
            break;
        case '6':
            u.ax = u.ax * cos(-1.0 * angle) + r.ax * sin(-1.0 * angle);
            u.ay = u.ay * cos(-1.0 * angle) + r.ay * sin(-1.0 * angle);
            u.az = u.az * cos(-1.0 * angle) + r.az * sin(-1.0 * angle);
            r = l.cross(u);
            break;

        default:
            break;
    }
}
void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
        case GLUT_KEY_UP:
            pos.x += l.ax * increment; //as l.x is negative incrementing value linearly will change position of eye
            pos.y += l.ay * increment; //as l.x is negative incrementing value linearly will change position of eye
            pos.z += l.az * increment;
            break;
        case GLUT_KEY_DOWN:
            pos.x -= l.ax * increment;
            pos.y -= l.ay * increment;
            pos.z -= l.az * increment;
            break;

        case GLUT_KEY_LEFT:
            pos.x -= r.ax * increment; //as right vector is perpendicular to up vector linear shifting with negative incrementing value linearly will change position of eye left side
            pos.y -= r.ay * increment;
            pos.z -= r.az * increment;
            break;
        case GLUT_KEY_RIGHT:
            pos.x += r.ax * increment; //as right vector is perpendicular to up vector linear shifting with positive incrementing value linearly will change position of eye right side
            pos.y += r.ay * increment;
            pos.z += r.az * increment;
            break;

        case GLUT_KEY_PAGE_UP:
            pos.x += u.ax * increment;
            pos.y += u.ay * increment;
            pos.z += u.az * increment; // shifting of position of eye on in z axis
            break;
        case GLUT_KEY_PAGE_DOWN:
            pos.x -= u.ax * increment;
            pos.y -= u.ay * increment;
            pos.z -= u.az * increment; // shifting of position of eye on in z axis
            break;

        default:
            break;
    }
}

int main(int argc, char **argv)
{
    int items ;
    string item_type;
    freopen("input.txt","r",stdin);

    cin>>recursionLevel ;
    cin>>width ;
    width = height;
    cin>>items ;

    for(int i=0;i<items;i++){
        cin>>item_type ;
        if(item_type=="pyramid"){
            Point low;
            cin>>low.x>>low.y>>low.z;
            double width, height;
            cin>>width>>height;
            Color color;
            cin>>color.red>>color.green>>color.blue;
            double am,df,spc, refl,shn;
            cin>>am>>df>>spc>>refl;
            cin>>shn;
            Pyramid pyramid(low,width,height,color);
            pyramids.push_back(pyramid);
        }

        if(item_type=="sphere"){
            Point center;
            cin>>center.x>>center.y>>center.z;
            double radius;
            cin>>radius;
            Color color;
            cin>>color.red>>color.green>>color.blue;
            double am,df,spc, refl,shn;
            cin>>am>>df>>spc>>refl;
            cin>>shn;
            Sphere sphere(center,radius,color);
            spheres.push_back(sphere);
        }
    }


    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutCreateWindow("Offline 3");
    init();
    glEnable(GL_DEPTH_TEST);
    glutDisplayFunc(display);
    glutIdleFunc(animate);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();
    return 0;
}
