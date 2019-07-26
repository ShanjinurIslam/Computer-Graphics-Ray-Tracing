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
#define pi acos(-1.0)

using namespace std ;

double plane_height=500;
double plane_width=500;
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

    Point operator-(Point v)
    {
        Point temp;
        temp.x = this->x - v.x;
        temp.y = this->y - v.y;
        temp.z = this->z - v.z;
        return temp;
    }

    Point &operator=(Point v)
    {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;

        return *this;
    }

    double distance(Point a){
        Point tem = *this - a ;
        return abs(sqrt(tem.x*tem.x+tem.y*tem.y+tem.z*tem.z)) ;
    }

    void print(){
        cout<<"( "<<x<<" , "<<y<<" , "<<z<<" )"<<endl;
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

Point lineParametric(Point p,Vector v,double t){
    Point new_point ;
    new_point.x = p.x + v.ax*t ;
    new_point.y = p.y + v.ay*t ;
    new_point.z = p.z + v.az*t ;

    return new_point ;
}

class Ray{
public:
    Point point ;
    Vector dir ;
    Ray(){
        this->point = Point() ;
        this->dir = Vector() ;
    }
    Ray(Point p,Vector v){
        this->point = p ;
        this->dir = v ;
    }

    Ray(Point a,Point b){
        dir = dir.generateVector(a,b) ;
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

    Color operator+(Color a){
        Color temp;
        temp.red = this->red+a.red;
        temp.green = this->green+a.green;
        temp.blue = this->blue+a.blue;
        return temp;
    }

    Color operator+(double a){
        Color temp;
        temp.red = this->red+a;
        temp.green = this->green+a;
        temp.blue = this->blue+a;
        return temp;
    }

    Color operator*(double a)
    {
        Color temp;
        temp.red = this->red*a;
        temp.green = this->green*a;
        temp.blue = this->blue*a;
        return temp;
    }

    Color &operator=(Color c)
    {
        this->red = c.red;
        this->green = c.green;
        this->blue = c.blue;

        return *this;
    }

    void printColor(){
        cout<<red*255<<" "<<green*255<<" "<<blue*255<<endl ;
    }
};

class Object{
public:
    int type ;
    Color color ;
    double a,d,s,r,shine;
    virtual double getT(Ray r) = 0;
    virtual double getColor(Ray r,Color *color,int level) = 0 ;
    virtual void draw() = 0 ;
    virtual Color getTileColor(Point p) = 0;
};

vector<Object*> objects ;
vector<Point> lights ;

class Triangle:public Object{
public:
    Point x,y,z ;
    Color color ;
    Triangle(){
    }

    Triangle(Point x,Point y,Point z,Color color,double a,double d,double s,double r,double shine){
        this->a = a ;
        this->d = d ;
        this->s = s ;
        this->r = r ;
        this->shine = shine ;
        this->type = 1 ;
        this->color = color ;
        this->x = x ;
        this->y = y ;
        this->z = z ;
    }

    Vector normal(){
        Vector n =  Vector().generateVector(x,y).cross(Vector().generateVector(x,z));
        n.normalize() ;
        return n ;
    }

    Vector normal(Vector r){
        Vector n =  Vector().generateVector(x,y).cross(Vector().generateVector(x,z));
        if(n.dot(r)>0) n = n*-1 ;
        n.normalize() ;
        return n ;
    }

    Vector getCross(Point a,Point b,Point c){
        Vector v1 = Vector().generateVector(a,b) ;
        Vector v2 = Vector().generateVector(a,c) ;

        Vector r = v1.cross(v2) ;

        return r ;
    }

    double getT(Ray r){
        Vector n = normal() ;
        double d = -(n.ax*x.x + n.ay*x.y + n.az*x.z) ;
        double t = -((d+n.dot(Vector().generateVector(Point(0,0,0),r.point)))/n.dot(r.dir));
        Point p = lineParametric(r.point,r.dir,t) ;

        Vector i,j,k ;

        i = getCross(x,y,p) ;
        j = getCross(y,z,p) ;
        k = getCross(z,x,p) ;

        if(i.dot(j)>=1 && i.dot(k)>=1) return t ;

        return -1 ;
    }

    double getColor(Ray ray,Color *out_color,int level){
        double t = getT(ray);
        if(t<=0) return -1;
        if(level==0)return t;
        *out_color = this->color * a ;
        Point ip = lineParametric(ray.point,ray.dir,t) ;
        Vector n = normal(ray.dir) ;

        for(int i=0;i<lights.size();i++){
            Vector l_dir = Vector().generateVector(ip,lights[i]);
            l_dir.normalize();
            Point s_p = lineParametric(ip,l_dir,1);
            Vector s_l = Vector().generateVector(s_p,lights[i]);
            s_l.normalize() ;
            Ray l_r(s_p,s_l) ;
            int touch = 0 ;
            double distance = s_p.distance(lights[i]);
            for(int k=0;k<objects.size();k++){
                double t = objects[k]->getT(l_r) ;
                if(t>0 && t<distance){
                    touch = 1 ;
                    break ;
                }
            }
            if(touch==1){
                Ray indt(lights[i],l_dir*-1);
                Vector light_ref = n*(2.0*indt.dir.dot(n)) - indt.dir ;
                double lambert = l_r.dir.dot(n);
                double phong = pow((ray.dir*-1).dot(light_ref),shine) ;

                lambert = max(lambert,0.0);
                phong = max(phong,0.0);

                *out_color = *out_color + (lambert*d + phong*s) ;
            }
        }

        if(level>0)
        {
            Vector ref =  ray.dir - n*(2.0*ray.dir.dot(n))  ;
            Point start = lineParametric(ip,ref,1) ;
            Point other = lineParametric(ip,ref,2);
            Vector v = Vector().generateVector(start,other) ;

            v.normalize() ;
            Ray reflectionRay(start,v);
            double t_min=INT_MAX;

            for(int k=0;k<objects.size();k++){
                Color *c = new Color();
                if(objects[k]->type==1){
                    double t = objects[k]->getColor(reflectionRay,c,recursionLevel-1);
                    if(t>0){
                        if(t_min>t){
                            t_min = t;
                            *out_color = *out_color + *c*r ;
                        }
                    }
                }
                if(objects[k]->type==2){
                    double t = objects[k]->getColor(reflectionRay,c,recursionLevel-1) ;
                    if(t>0){
                        if(t_min>t){
                            t_min = t;
                            Point p = lineParametric(reflectionRay.point,reflectionRay.dir,t);
                            if((p.x>=-750 && p.x<=750) && (p.y>=-750 && p.y<=750)){
                                *out_color = *out_color + *c*r ;
                            }
                        }
                    }
                }
                delete c;
            }

        }
        return t ;
    }

    void draw(){
        glColor3f(color.red,color.green,color.blue);
        glBegin(GL_TRIANGLES);
        glVertex3d(x.x,x.y,x.z);
        glVertex3d(y.x,y.y,y.z);
        glVertex3d(z.x,z.y,z.z);
        glEnd();
    }

    Color getTileColor(Point p){
        return Color(0,0,0);
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

    ~Square(){
        delete[] corner;
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

class Sphere:public Object{
public:
    Point center ;
    double radius;

    Sphere(Point center,double radius,Color color,double a,double d,double s,double ref,double shine){
        this->a = a ;
        this->d = d ;
        this->s = s ;
        this->r = ref ;
        this->shine = shine ;
        this->center = center ;
        this->radius = radius ;
        this->color = color ;
        this->type = 1 ;
    }

    double getT(Ray r){

        Vector r0 = Vector().generateVector(center,r.point);
        Vector rd = r.dir ;
        double a = 1 ;
        double b = 2*r0.dot(rd);
        double c = r0.dot(r0)- radius*radius;
        double d = b*b - (4*a*c) ;

        if(d<0){
            return -1 ;
        }

        double alpha = (-b + sqrt(d))/2 ;
        double beta = (-b - sqrt(d))/2 ;

        double t = min(alpha,beta);

        if(t<=0){
            t = max(alpha,beta);
        }

        return t ;
    }

    double getColor(Ray ray,Color *out_color,int level){
        double t = getT(ray);
        *out_color = this->color * a ;
        if(t<=0) return -1;

        if(level==0)return t;

        Point ip = lineParametric(ray.point,ray.dir,t) ;
        Vector n = Vector().generateVector(center,ip);
        n.normalize() ;

        for(int i=0;i<lights.size();i++){
            Vector l_dir = Vector().generateVector(ip,lights[i]);
            l_dir.normalize();
            Point s_p = lineParametric(ip,l_dir,1);
            Vector s_l = Vector().generateVector(s_p,lights[i]);
            s_l.normalize() ;
            Ray l_r(s_p,s_l) ;
            int touch = 0 ;
            double distance = s_p.distance(lights[i]);
            for(int k=0;k<objects.size();k++){
                double t = objects[k]->getT(l_r) ;
                if(t>0 && t<distance){
                    touch = 1 ;
                    break ;
                }
            }
            if(touch==1){
                Ray indt(lights[i],l_dir*-1);
                Vector light_ref = n*(2.0*indt.dir.dot(n)) - indt.dir ;
                double lambert = l_r.dir.dot(n);
                double phong = pow((ray.dir*-1).dot(light_ref),shine) ;

                lambert = max(lambert,0.0);
                phong = max(phong,0.0);

                *out_color = *out_color + (lambert*d + phong*s) ;
            }
        }

        if(level>0)
        {
            Vector ref = ray.dir - n*(2.0*ray.dir.dot(n))   ;
            Point start = lineParametric(ip,ref,1) ;
            Point other = lineParametric(ip,ref,2);
            Vector v = Vector().generateVector(start,other) ;
            v.normalize() ;
            Ray reflectionRay(start,v);
            double t_min=INT_MAX;

            for(int k=0;k<objects.size();k++){
                Color *c = new Color();
                if(objects[k]->type==1){
                    double t = objects[k]->getColor(reflectionRay,c,recursionLevel-1);
                    if(t>0){
                        if(t_min>t){
                            t_min = t;
                            *out_color = *out_color + *c*r ;
                        }
                    }
                }
                if(objects[k]->type==2){
                    double t = objects[k]->getColor(reflectionRay,c,recursionLevel-1) ;
                    if(t>0){
                        if(t_min>t){
                            t_min = t;
                            Point p = lineParametric(reflectionRay.point,reflectionRay.dir,t);
                            if((p.x>=-750 && p.x<=750) && (p.y>=-750 && p.y<=750)){
                                *out_color = *out_color + *c*r ;
                            }
                        }
                    }
                }
                delete c;
            }

        }
        return t ;
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

    Color getTileColor(Point p){
        return Color(0,0,0);
    }
};

class ChessBoard:public Object{
private:
    int color_map[51][51];
public:
    ChessBoard(){
        this->type = 2 ;
        int flag = 1 ;
        for(int i=0;i<=50;i++){
            for(int j=0;j<=50;j++){
                color_map[i][j] = flag%2 ;
                flag++;
            }
        }
    }

    double getT(Ray r){
        return -(Vector(0,0,1).dot(Vector().generateVector(Point(0,0,0),r.point))/Vector(0,0,1).dot(r.dir));
    }

    double getColor(Ray ray,Color *out_color,int level){
        double t = getT(ray);

        if(t<=0) return -1;
        if(level==0)return t;

        Point ip = lineParametric(ray.point,ray.dir,t) ;

        Vector n = Vector(0,0,1);

        if((ip.x>=-750 && ip.x<=750) && (ip.y>=-750 && ip.y<=750)){
            *out_color = getTileColor(ip) * 0.4 ;
            for(int i=0;i<lights.size();i++){
                Vector l_dir = Vector().generateVector(ip,lights[i]);
                l_dir.normalize();
                Point s_p = lineParametric(ip,l_dir,1);
                Vector s_l = Vector().generateVector(s_p,lights[i]);
                s_l.normalize() ;
                Ray l_r(s_p,s_l) ;
                int touch = 0 ;
                double distance = s_p.distance(lights[i]);
                for(int k=0;k<objects.size();k++){
                    double t = objects[k]->getT(l_r) ;
                    if(t>0 && t<distance){
                        touch = 1 ;
                        break ;
                    }
                }
                if(touch==1){
                    Ray indt(lights[i],l_dir*-1);
                    Vector light_ref =  n*(2.0*indt.dir.dot(n)) - indt.dir ;
                    double lambert = l_r.dir.dot(n);
                    double phong = pow((ray.dir*-1).dot(light_ref),5) ;

                    lambert = max(lambert,0.0);
                    phong = max(phong,0.0);

                    *out_color = *out_color + (lambert*0.2 + phong*0.2) ;
                }
            }
        }



        if(level>0)
        {
            Vector ref = n*(2.0*ray.dir.dot(n)) - ray.dir ;
            Point start = lineParametric(ip,ref,1) ;
            Point other = lineParametric(ip,ref,2);
            Vector v = Vector().generateVector(start,other) ;
            v.normalize() ;
            Ray reflectionRay(start,v);
            double t_min=INT_MAX;

            for(int k=0;k<objects.size();k++){
                Color *c = new Color();
                if(objects[k]->type==1){
                    double t = objects[k]->getColor(reflectionRay,c,recursionLevel-1);
                    if(t>0){
                        if(t_min>t){
                            t_min = t;
                            *out_color = *out_color + *c*r ;
                        }
                    }
                }
                if(objects[k]->type==2){
                    double t = objects[k]->getColor(reflectionRay,c,recursionLevel-1) ;
                    if(t>0){
                        if(t_min>t){
                            t_min = t;
                            Point p = lineParametric(reflectionRay.point,reflectionRay.dir,t);
                            if((p.x>=-750 && p.x<=750) && (p.y>=-750 && p.y<=750)){
                                *out_color = *out_color + *c*r ;
                            }
                        }
                    }
                }
                delete c;
            }
        }

        return t ;
    }



    Color getTileColor(Point p){
        if(color_map[(int)floor(p.x/30.0)+25][(int)floor(p.y/30.0)+25]){
            return Color(0,0,0) ;
        }
        else{
            return Color(1,1,1) ;
        }
    }

    void draw(){
        int flag=1;
        for(double i=-30*25;i<=30*25;i+=30){
            for(double j=-30*25;j<=30*25;j+=30){
                if(color_map[(int)(i/30.0)+25][(int)(j/30.0)+25]){
                    Square(Point(i,j,0),30,Color(0,0,0)).draw();
                }
                else{
                    Square(Point(i,j,0),30,Color(1,1,1)).draw();
                }
                flag++;
            }
        }
    }
};

ChessBoard *chessboard ;

void init()
{
    cameraAngle=80;
    pos = Point(50, 50, 50);
    u = Vector(0, 0, 1); //z axis is up vector
    l = Vector(-1 / (sqrt(2)), -1 / sqrt(2), 0);
    r = Vector(-1 / (sqrt(2)), 1 / sqrt(2), 0);
    angle = acos(-1.0) / 30; //3 degree angle change
    increment = 5 ;
    chessboard = new ChessBoard();
    objects.push_back(chessboard);
    glClearColor(0, 0, 0, 0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(cameraAngle, 1, 1, 1000.0);
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

    for(int i=0;i<objects.size();i++){
        objects[i]->draw() ;
    }

    for(int i=0;i<lights.size();i++){
        glColor3f(1.0,1.0,1.0);
        glPointSize(2);
        glBegin(GL_POINTS);
        glVertex3d(lights[i].x,lights[i].y,lights[i].z);
        glEnd();
    }

    glPopMatrix();

    glutSwapBuffers();
}
void animate()
{
    glutPostRedisplay();
}

void generateRayTracedImage(){
    double near_plane_distance = (plane_height/2.0)*0.8391;
    Vector new_l = l*near_plane_distance ;
    Vector new_r = r*(-1.0*(plane_height/2.0)) ;
    Vector new_u = u*(plane_width/2.0) ;
    Vector resultant = new_l+new_r+new_u ;
    Point top_left = lineParametric(pos,resultant,1);

    double pixel_x = plane_height/height ;
    double pixel_y = plane_width/width ;

    int image_height = (int) height ;
    int image_width = (int) width ;

    bitmap_image image(image_height,image_width);

    for(int i=0;i<image_height;i++){
        for(int j=0;j<image_width;j++){
            Point point ;
            Vector n_u = u*(-1.0*i*pixel_x);
            Vector n_r = r*(j*pixel_y);
            Vector n_resultant = n_u + n_r ;
            point = lineParametric(top_left,n_resultant,1);
            Vector v = Vector().generateVector(pos,point);
            v.normalize();
            Ray ray(pos,v);
            double t_min = INT_MAX ;

            for(int k=0;k<objects.size();k++){
                Color *c = new Color();
                if(objects[k]->type==1){
                    double t = objects[k]->getColor(ray,c,recursionLevel);
                    if(t>0){
                        if(t_min>t){
                            t_min = t;
                            image.set_pixel(j, i,c->red*255,c->green*255,c->blue*255);
                        }
                    }
                }
                if(objects[k]->type==2){
                    double t = objects[k]->getColor(ray,c,recursionLevel) ;
                    if(t>0){
                        if(t_min>t){
                            t_min = t;
                            Point p = lineParametric(ray.point,ray.dir,t);
                            if((p.x>=-750 && p.x<=750) && (p.y>=-750 && p.y<=750)){
                                image.set_pixel(j, i, c->red*255,c->green*255,c->blue*255);
                            }
                        }
                    }
                }
                delete c;
            }

            /*for(int k=0;k<spheres.size();k++){
             double t = spheres[k].getColor(ray);
             if(t>0){
             if(t_min>t){
             t_min = t;
             Color c = spheres[k].color ;
             image.set_pixel(j, i, (int) c.red*255, (int) c.green*255,(int) c.blue*255);
             }
             }
             }


             for(int k=0;k<pyramids.size();k++){
             double t = pyramids[k].getColor(ray);
             if(t>0){
             if(t_min>t){
             t_min = t;
             Color c = pyramids[k].color ;
             image.set_pixel(j, i, (int) c.red*255, (int) c.green*255,(int) c.blue*255);
             }
             }
             }

             double t = chessboard.getColor(ray) ;
             if(t>0){
             if(t_min>t){
             t_min = t;
             Color c = chessboard.getTileColor(lineParametric(ray.point,ray.dir,t));
             image.set_pixel(j, i, (int) c.red*255, (int) c.green*255,(int) c.blue*255);
             }
             }*/
        }
    }

    image.save_image("out.bmp");
    printf("Image generated\n");
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
        case '0':
            generateRayTracedImage() ;
            break ;
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
    height = width;
    cin>>items ;

    for(int i=0;i<items;i++){
        double a,d,s,r,shine;
        cin>>item_type ;
        if(item_type=="pyramid"){
            Point low;
            cin>>low.x>>low.y>>low.z;
            double width, height;
            cin>>width>>height;
            Color color;
            cin>>color.red>>color.green>>color.blue;
            cin>>a>>d>>s>>r;
            cin>>shine;
            Square bottom ;
            Point top = low + Point(width/2,width/2,height);
            bottom = Square(low,width,color) ;

            Point *corners = new Point[4] ;
            corners = bottom.getCornerPoints() ;

            for(int i=0;i<4;i++){
                Triangle *t = new Triangle(corners[i],corners[(i+1)%4],top,color,a,d,s,r,shine) ;
                objects.push_back(t);
            }
        }

        if(item_type=="sphere"){
            Point center;
            cin>>center.x>>center.y>>center.z;
            double radius;
            cin>>radius;
            Color color;
            cin>>color.red>>color.green>>color.blue;
            cin>>a>>d>>s>>r;
            cin>>shine;
            Sphere *sphere = new Sphere(center,radius,color,a,d,s,r,shine);
            objects.push_back(sphere);
        }
    }

    cin>>items ;
    for(int i=0;i<items;i++){
        Point low;
        cin>>low.x>>low.y>>low.z;
        lights.push_back(Point(low.x,low.y,low.z));
        lights[i].print() ;
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

