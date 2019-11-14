#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    float sum=0;
      for(int w=0;w<im.w;w++){
        for(int h=0;h<im.h;h++){
            for(int c=0;c<im.c;c++){
                sum+=get_pixel(im,w,h,c);
            }
        }
    }
    for(int w=0;w<im.w;w++){
        for(int h=0;h<im.h;h++){
            for(int c=0;c<im.c;c++){
                set_pixel(im,w,h,c,get_pixel(im,w,h,c)/sum);

        }
    }
    }
}

image make_box_filter(int w)
{
    // TODO
    image re=make_image(w,w,1);
    for (int x=0;x<w;x++){
        for(int y=0;y<w;y++){
            set_pixel(re,x,y,0,1.0);
        }
    }
    l1_normalize(re);
    return re;
}

image convolve_image(image im, image filter, int preserve)
{   assert(im.c == filter.c || filter.c == 1);
    image re;

    if(preserve==1){
        re=make_image(im.w,im.h,im.c);
    }
    else{
        re=make_image(im.w,im.h,1);
    }
    for (int x=0;x<im.w;x++){
            for(int y=0; y<im.h;y++){
                for(int c=0;c<im.c;c++){
                    float v=0;
                    for(int fx=0;fx<filter.w;fx++){
                        for(int fy=0;fy<filter.h;fy++){
                            v+=get_pixel(im,x+fx-filter.w/2,y+fy-filter.h/2,c)*get_pixel(filter,fx,fy,0);
                        }
                    }
                     if(preserve==1){
                        set_pixel(re,x,y,c,v);
                     }
                    else{
                        float v0=get_pixel(re,x,y,0)+v;
                        set_pixel(re,x,y,0,v0);
                    }
                }
        }
    }
    return re;
}

image make_highpass_filter()
{
    // TODO
    image re=make_box_filter(3);
    re.data[0]=0;
    re.data[1]=-1;
    re.data[2]=0;
    re.data[3]=-1;
    re.data[4]=4;
    re.data[5]=-1;
    re.data[6]=0;
    re.data[7]=-1;
    re.data[8]=0;

    return re;
}

image make_sharpen_filter()
{
    image re=make_box_filter(3);
    re.data[0]=0;
    re.data[1]=-1;
    re.data[2]=0;
    re.data[3]=-1;
    re.data[4]=5;
    re.data[5]=-1;
    re.data[6]=0;
    re.data[7]=-1;
    re.data[8]=0;

    return re;
}

image make_emboss_filter()
{
    image re=make_box_filter(3);
    re.data[0]=-2;
    re.data[1]=-1;
    re.data[2]=0;
    re.data[3]=-1;
    re.data[4]=1;
    re.data[5]=1;
    re.data[6]=0;
    re.data[7]=1;
    re.data[8]=2;

    return re;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: When we do emboss and sharpen, we expect the output stand in original color. Therefore, sharp and emboss filter need preserving.
// On the other hand, in highpass filter, we do not preserve the color, so highpass does not need preserving.

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Yes, we have to do post-processing on each filter. If we do not do post-processing, the picture will be out of constraints.
// As result, the outcome will perform weird by out of constraints.

image make_gaussian_filter(float sigma)
{
    // TODO
    //  the kernel be 6 times the size of sigma.odd number, so make it be the next highest odd integer from 6x sigma.
    int s=(int) roundf(sigma*6)+1;
    int w=0;
    if(s%2!=0){
        w=s;
    }
    else{
        w=s+1;
    }
    image re=make_image(w,w,1);
    int center=w/2;
    for(int x=0;x<re.w;x++){
        float v=0;
        for(int y=0;y<re.w;y++)
            {
                int y1=y-center;
                int x1=x-center;
                v=1/(TWOPI*sigma*sigma)*exp(-1*(x1*x1+y1*y1)/(2*sigma*sigma));
                set_pixel(re,x,y,0,v);
            }

    }
    l1_normalize(re);

    return re;
}

image add_image(image a, image b)
{
    // TODO
     assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image re=make_image(a.w,a.h,a.c);
    for (int x=0;x<a.w;x++){
        for(int y=0; y<a.w;y++){
            for(int c=0;c<a.c;c++){
                set_pixel(re,x,y,c,get_pixel(a,x,y,c)+get_pixel(b,x,y,c));
            }
        }
    }
    return re;
}

image sub_image(image a, image b)
{
    // TODO
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image re=make_image(a.w,a.h,a.c);
    for (int x=0;x<a.w;x++){
        for(int y=0; y<a.w;y++){
            for(int c=0;c<a.c;c++){
                set_pixel(re,x,y,c,get_pixel(a,x,y,c)-get_pixel(b,x,y,c));
            }
        }
    }
    return re;
}

image make_gx_filter()
{
  image re = make_box_filter(3);
  re.data[0] = -1;
  re.data[1] = 0;
  re.data[2] = 1;
  re.data[3] = -2;
  re.data[4] = 0;
  re.data[5] = 2;
  re.data[6] = -1;
  re.data[7] = 0;
  re.data[8] = 1;
  return re;
}

image make_gy_filter()
{
image re = make_box_filter(3);
  re.data[0] = -1;
  re.data[1] = -2;
  re.data[2] = -1;
  re.data[3] = 0;
  re.data[4] = 0;
  re.data[5] = 0;
  re.data[6] = 1;
  re.data[7] = 2;
  re.data[8] = 1;
  return re;
}

void feature_normalize(image im)
{
    // TODO
    float minn=999999;
    float Maxx=-1;
    for(int x=0;x<im.w;x++){
        for(int y=0;y<im.h;y++){
            for(int c=0;c<im.c;c++){
                if(get_pixel(im,x,y,c)>Maxx){
                    Maxx=get_pixel(im,x,y,c);
                }
                 if(get_pixel(im,x,y,c)<minn){
                    minn=get_pixel(im,x,y,c);
                }
            }
        }
    }
    if(Maxx-minn!=0){
        for(int x=0;x<im.w;x++){
            for(int y=0;y<im.h;y++){
                for(int c=0;c<im.c;c++){
                    float range=Maxx-minn;
                    set_pixel(im,x,y,c,get_pixel(im,x,y,c)/range);
                    }
                }
            }
        }
else{
    for(int x=0;x<im.w;x++){
        for(int y=0;y<im.h;y++){
            for(int c=0;c<im.c;c++){
                 set_pixel(im,x,y,c,0);
            }
        }
    }
}

}

image *sobel_image(image im)
{
    // TODO
    image fgx=make_gx_filter();
    image fgy=make_gy_filter();
    image gx=convolve_image(im,fgx,0);
    image gy=convolve_image(im,fgy,0);
    image *re=calloc(2, sizeof(image));

    re[0]=make_image(im.w,im.h,1);
    re[1]=make_image(im.w,im.h,1);
    for(int x=0;x<im.w;x++){
        for(int y=0;y<im.h;y++){
            float magn=sqrt(pow(get_pixel(gx,x,y,0),2)+pow(get_pixel(gy,x,y,0),2));
            float gradient=atan2(get_pixel(gy,x,y,0),get_pixel(gx,x,y,0));
            set_pixel(re[0],x,y,0,magn);
            set_pixel(re[1],x,y,0,gradient);
        }
    }
    return re;
}

image colorize_sobel(image im)
{
    image *s = sobel_image(im);
    feature_normalize(s[0]);
    feature_normalize(s[1]);
    image re = make_image(im.w, im.h, 3);
    for(int x=0;x<im.w;x++){
        for(int y=0;y<im.h;y++){
            set_pixel(re,x,y,0,get_pixel(s[1],x,y,0));
            set_pixel(re,x,y,1,get_pixel(s[0],x,y,0));
            set_pixel(re,x,y,2,get_pixel(s[0],x,y,0));
        }
    }
    hsv_to_rgb(re);
    return  re;
}
