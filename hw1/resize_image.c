#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
    return get_pixel(im,(int)(x+0.5),(int)(y+0.5),c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)

    image nimage=make_image(w,h,im.c);
    float s_w=(float)im.w/w;
    float s_h=(float)im.h/h;
    for(int i=0;i<w;i++){
        for(int j=0;j<h;j++){
            for(int c=0;c<im.c;c++){
                  float o_x=s_w*i-0.5+0.5*s_w;
                  float o_y=s_h*j-0.5+0.5*s_h;
                  set_pixel(nimage,i,j,c,nn_interpolate(im,o_x,o_y,c));
            }
            }
    }
    return nimage;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
            int o_w=(int)x;
            int o_h=(int)y;
            int o_w_1=(int)(x+1);
            int o_h_1=(int)(y+1);
            float t_w=(float)x;
            float t_h=(float)y;
            if(o_w_1==im.w){
                o_w_1=o_w;
            }
            if(o_h_1==im.h){
                o_h_1=o_h;
            }
            //Area total

            //each weight
            float a1=fabs((o_w_1-t_w)* (o_h_1-t_h));
            float a2=fabs((o_w-t_w)* (o_h_1-t_h));
            float a3=fabs((o_w_1-t_w)* (o_h-t_h));
            float a4=fabs((o_w-t_w)* (o_h-t_h));
            float area=a1+a2+a3+a4;
            //each value
            float v1=get_pixel(im,o_w,o_h,c);
            float v2=get_pixel(im,o_w_1,o_h,c);
            float v3=get_pixel(im,o_w,o_h_1,c);
            float v4=get_pixel(im,o_w_1,o_h_1,c);
            float value=(v1*a1+v2*a2+v3*a3+v4*a4)/area;
    return value;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    image nimage=make_image(w,h,im.c);
    float s_w=(float)im.w/w;
    float s_h=(float)im.h/h;
    for(int i=0;i<w;i++){
        for(int j=0;j<h;j++){
            for(int c=0;c<im.c;c++){
            float o_x=s_w*i-0.5+0.5*s_w;
            float o_y=s_h*j-0.5+0.5*s_h;
            set_pixel(nimage,i,j,c,bilinear_interpolate(im,o_x,o_y,c));
            }
        }
    }
    return nimage;
}


