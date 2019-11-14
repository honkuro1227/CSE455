#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    if (x < 0) {
	x = 0;
    }
    if (x >= im.w) {
	x = im.w - 1;
    }
    if (y < 0) {
	y = 0;
    }
    if (y >= im.h) {
	y = im.h - 1;
    }
    if (c < 0) {
	c = 0;
    }
    if (c >= im.c) {
	c = im.c - 1;
    }
    return im.data[x + y*im.w + c*im.w*im.h];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if ((x >= 0 && x < im.w) && (y >= 0 && y < im.h) && (c >= 0 && c < im.c)) {
	im.data[x + y*im.w + c*im.w*im.h] = v;
    }
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, im.w*im.h*im.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    for (int x = 0; x < im.w; x++) {
	for (int y = 0; y < im.h; y++) {
	    gray.data[x + y*im.w] = 0.299*im.data[x + y*im.w] + 0.587*im.data[x + y*im.w + im.w*im.h] + 0.114*im.data[x + y*im.w + 2*im.w*im.h];
	}
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    if (c >= 0 && c < im.c) {
	for (int x = 0; x < im.w; x++) {
	    for (int y = 0; y < im.h; y++) {
		im.data[x + y*im.w + c*im.w*im.h] += v;
	    }
	}
    }
}

void clamp_image(image im)
{
    for (int i = 0; i < im.w*im.h*im.c; i++) {
	if (im.data[i] < 0) {
	    im.data[i] = 0;
	}
	if (im.data[i] > 1) {
	    im.data[i] = 1;
	}
    }
}

// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    float r, g, b, min, max, c, h, s, v;
    for (int x = 0; x < im.w; x++) {
	for (int y = 0; y < im.h; y++) {
	    r = im.data[x + y*im.w];
	    g = im.data[x + y*im.w + im.w*im.h];
	    b = im.data[x + y*im.w + 2*im.w*im.h];
	    min = three_way_min(r, g, b);
	    max = three_way_max(r, g, b);
	    c = max - min;
	    v = max;
	    if (v != 0) {
		s = c / v;
	    }
	    else {
	    	s = 0;
	    }
	    if (c != 0) {
	    	if (v == r) {
	    	    h = (g - b) / c;
	    	}
	    	else if (v == g) {
	    	    h = ((b - r) / c) + 2;
	    	}
	    	else {
	    	    h = ((r - g) / c) + 4;
	    	}
	    }
	    else {
	    	h = 0;
	    }
	    h /= 6;
	    if (h < 0) {
	    	h++;
	    }
	    im.data[x + y*im.w] = h;
	    im.data[x + y*im.w + im.w*im.h] = s;
	    im.data[x + y*im.w + 2*im.w*im.h] = v;
	}
    }
}

void hsv_to_rgb(image im)
{
    float h, s, v, c, max, min, r, g, b, h_temp;
    for (int x = 0; x < im.w; x++) {
	for (int y = 0; y < im.h; y++) {
	    h = im.data[x + y*im.w];
	    s = im.data[x + y*im.w + im.w*im.h];
	    v = im.data[x + y*im.w + 2*im.w*im.h];
	    c = s * v;
	    max = v;
	    if (v != c) {
		min = v - c;
	    }
	    else {
	    	min = 0.0;
	    }
	    h_temp = h * 6;
	    if (c == 0) {
	    	r = v;
	    	g = v;
	    	b = v;
	    }
	    else if (h_temp > 5 && h_temp < 6) {
	    	r = max;
	    	g = min;
	    	b = ((((h_temp /  6) - 1) * 6 * c) - g) * -1;
	    }
	    else if (h_temp == 5) {
	    	r = max;
	    	g = min;
	    	b = max;
	    }
            else if (h_temp < 5 && h_temp > 4) {
	    	g = min;
	    	r = (h_temp - 4) * c + g;
	    	b = max;
	    }
	    else if (h_temp == 4) {
	    	r = min;
	    	g = min;
	    	b = max;
	    }
	    else if (h_temp < 4 && h_temp > 3) {
	    	r = min;
	    	g = (((h_temp - 4) * c) - r) * -1;
	    	b = max;
	    }
            else if (h_temp == 3) {
	    	r = min;
	    	g = max;
	    	b = max;
	    }
            else if (h_temp < 3 && h_temp > 2) {
	    	r = min;
	    	g = max;
	    	b = ((h_temp - 2) * c) + r;
	    }
	    else if (h_temp == 2) {
	    	r = min;
	    	g = max;
	    	b = min;
	    }
            else if (h_temp < 2 && h_temp > 1) {
	    	g = max;
	    	b = min;
	    	r = (((h_temp - 2) * c) - b) * -1;
	    }
            else if (h_temp == 1) {
	    	r = max;
	    	g = max;
	    	b = min;
	    }
             else if (h_temp < 1 && h_temp > 0) {
	    	r = max;
	    	b = min;
	    	g = (h_temp * c) + b;
	    }
            else {
	    	r = max;
	    	g = min;
	    	b = min;
	    }
	    im.data[x + y*im.w] = r;
	    im.data[x + y*im.w + im.w*im.h] = g;
	    im.data[x + y*im.w + 2*im.w*im.h] = b;
	}
    }
}

void scale_image(image im, int c, float v)
{
    if (c >= 0 && c < im.c) {
	for (int x = 0; x < im.w; x++) {
	    for (int y = 0; y < im.h; y++) {
		im.data[x + y*im.w + c*im.w*im.h] *= v;
	    }
	}
    }
}
//0.8  RGB_to_HCL
void RGB_to_HCL(image im){
    //Assume D65 in obverse 2 degree, and Linear RGB
     float D65_X=0.95047;
     float D65_Y=1.00000;
     float D65_Z=1.08883;
     for (int col=0;col<im.w;col++){
        for(int row=0;row<im.h;row++){
            float R=get_pixel(im,col,row,0);
            float G=get_pixel(im,col,row,1);
            float B=get_pixel(im,col,row,2);
            //RGBtoXYZ
            float X = 0.4124564*R + 0.3575761*G + 0.1804375*B;
            float Y = 0.2126729*R + 0.7151522*G + 0.0721750*B;
            float Z = 0.0193339*R + 0.1191920*G + 0.9503041*B;
            //XYZtoCIELUV

            float L;
            float u;
            float v;

            if((Y/D65_Y)<=pow(6.0/29.0,3))
            {
                L=pow((29.0/3.0),3)*(Y/D65_Y);
            }
            if((Y/D65_Y)>pow(6.0/29.0,3)){
                L=1.16*pow((Y/D65_Y),1.0/3.0)-0.16;
            }

            float ua=4.0*X/(X+15.0*Y+3.0*Z);
            float va=9.0*Y/(X+15.0*Y+3.0*Z);
            float un=4.0*D65_X/(D65_X+15.0*D65_Y+3.0*D65_Z);
            float vn=9.0*D65_Y/(D65_X+15.0*D65_Y+3.0*D65_Z);
             u=13.0*L*(ua-un);
             v=13.0*L*(va-vn);
            //CIELUVtoCIELCH

            float H=atan2(v,u);
            if(H<0){
                H=H/6.0+1;
            }
            else{
                H=H/6.0;
            }
            float C=pow((pow(v,2)+pow(u,2)),0.5);
            set_pixel(im,col,row,0,L);
            set_pixel(im,col,row,1,C);
            set_pixel(im,col,row,2,H);
}
}
}
//0.8HCL_to_RGB
void HCL_to_RGB(image im){
    float D65_X=0.95047;
    float D65_Y=1.00000;
    float D65_Z=1.08883;
    float un=4.0*D65_X/(D65_X+15.0*D65_Y+3.0*D65_Z);
    float vn=9.0*D65_Y/(D65_X+15.0*D65_Y+3.0*D65_Z);
    for (int col=0;col<im.w;col++)
     {
        for(int row=0;row<im.h;row++)
        {   //CIELCHtoCIELuv
            float L=get_pixel(im,col,row,0);
            float C=get_pixel(im,col,row,1);
            float H=get_pixel(im,col,row,2);
            float T=6.0*H;

            float u=C*cos(T);
            float v=C*sin(T);

            //CIELUVtoXYZ
            float X,Y,Z;

            if(L<=0.08){
                Y=D65_Y*L*100.0*pow(3.0/29.0,3);
            }
            if(L>0.08){
                Y=D65_Y*pow((L+0.16)/1.16,3);
            }
            if(L!=0.0){
            float ua=u/(13.0*L)+un;
            float va=v/(13.0*L)+vn;
            X=Y*9.0*ua/(4.0*va);
            Z=Y*(12.0-3.0*ua-20.0*va)/(4.0*va);
            }
            else{
            X=0.0;
            Y=0.0;
            }
            //CIEXYZtoRGB
            float R = 3.2404542*X - 1.5371385*Y - 0.4985314*Z;
            float G = -0.9692660*X + 1.8760108*Y + 0.0415560*Z;
            float B = 0.0556434*X - 0.2040259*Y + 1.0572252*Z;
            set_pixel(im,col,row,0,R);
            set_pixel(im,col,row,1,G);
            set_pixel(im,col,row,2,B);
        }
     }
}



