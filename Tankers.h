
#include <stdio.h>
#include <stdlib.h> 
#include <strings.h>
#include <math.h> 
#include "lodepng.h" 

// принимаем на вход: имя файла, указатели на int для хранения прочитанной ширины и высоты картинки
// возвращаем указатель на выделенную память для хранения картинки
// Если память выделить не смогли, отдаем нулевой указатель и пишем сообщение об ошибке
unsigned char* load_png(const char* filename, unsigned int* width, unsigned int* height) 
{
  unsigned char* image = NULL; 
  int error = lodepng_decode32_file(&image, width, height, filename);
  if(error != 0) {
    printf("error %u: %s\n", error, lodepng_error_text(error)); 
  }
  return (image);
}

// принимаем на вход: имя файла для записи, указатель на массив пикселей,  ширину и высоту картинки
// Если преобразовать массив в картинку или сохранить не смогли,  пишем сообщение об ошибке
void write_png(const char* filename, const unsigned char* image, unsigned width, unsigned height)
{
  unsigned char* png;
  long unsigned int pngsize;
  int error = lodepng_encode32(&png, &pngsize, image, width, height);
  if(error == 0) {
      lodepng_save_file(png, pngsize, filename);
  } else { 
    printf("error %u: %s\n", error, lodepng_error_text(error));
  }
  free(png);
}


// вариант огрубления серого цвета в ЧБ 
void contrast(unsigned char *col, int bw_size)
{ 
    int i; 
    for(i=0; i < bw_size; i++)
    {
        if(col[i] < 55)
        col[i] = 0; 
        if(col[i] > 195)
        col[i] = 255;
    } 
    return; 
} 

// Гауссово размыттие
void Gauss_blur(unsigned char *col, unsigned char *blr_pic, int width, int height)
{ 
    int i, j; 
    for(i=1; i < height-1; i++) 
        for(j=1; j < width-1; j++)
        { 
            blr_pic[width*i+j] = 0.084*col[width*i+j] + 0.084*col[width*(i+1)+j] + 0.084*col[width*(i-1)+j]; 
            blr_pic[width*i+j] = blr_pic[width*i+j] + 0.084*col[width*i+(j+1)] + 0.084*col[width*i+(j-1)]; 
            blr_pic[width*i+j] = blr_pic[width*i+j] + 0.063*col[width*(i+1)+(j+1)] + 0.063*col[width*(i+1)+(j-1)]; 
            blr_pic[width*i+j] = blr_pic[width*i+j] + 0.063*col[width*(i-1)+(j+1)] + 0.063*col[width*(i-1)+(j-1)]; 
        } 
   return; 
} 

//  Место для экспериментов
void color(unsigned char *blr_pic, unsigned char *res, int size)
{ 
  int i;
    for(i=1;i<size;i++) 
    { 
        res[i*4]=40+blr_pic[i]+0.35*blr_pic[i-1]; 
        res[i*4+1]=65+blr_pic[i]; 
        res[i*4+2]=170+blr_pic[i]; 
        res[i*4+3]=255; 
    } 
    return; 
}

int count_tankers(unsigned char* bw_pic, int startx, int starty, int width, int height, int realwidth) {
    int count = 0;
    for(int i = starty; i <starty + height; i++) {
        for(int j = startx; j <startx + width; j++) {
            int index = i * realwidth + j;
            if(bw_pic[index] == 255) {
                count++;
            }
        }
    }
    return count;
}
  
int main() 
{ 
    const char* filename = "skull.png"; 
    unsigned int width, height;
    int size;
    int bw_size;

    int all_tankers = 0;
    
    // Прочитали картинку
    unsigned char* picture = load_png("skull.png", &width, &height); 
    if (picture == NULL)
    { 
        printf("Problem reading picture from the file %s. Error.\n", filename); 
        return -1; 
    }

    size = width * height * 4;
    bw_size = width * height;
    
    unsigned char* bw_pic = (unsigned char*)malloc(bw_size*sizeof(unsigned char)); 
    unsigned char* blr_pic = (unsigned char*)malloc(bw_size*sizeof(unsigned char)); 
    unsigned char* finish = (unsigned char*)malloc(size*sizeof(unsigned char)); 

    for(int i = 0; i < bw_size; i++) {
        unsigned char r = picture[i*4];
        unsigned char g = picture[i*4+1];
        unsigned char b = picture[i*4+2];
        bw_pic[i] = 0.299 * r + 0.587 * g + 0.114 * b;
    }

    // Например, поиграли с  контрастом
    contrast(bw_pic, bw_size); 

    all_tankers += count_tankers(bw_pic, 400, 3, 135, 181, width);
    all_tankers += count_tankers(bw_pic, 379, 188, 316, 49, width);
    all_tankers += count_tankers(bw_pic, 411, 237, 15, 7, width);
    all_tankers += count_tankers(bw_pic, 380, 245, 463, 75, width);
    all_tankers += count_tankers(bw_pic, 396, 320, 462, 78, width);
    all_tankers += count_tankers(bw_pic, 400, 401, 374, 21, width);
    all_tankers += count_tankers(bw_pic, 425, 423, 399, 29, width);
    all_tankers += count_tankers(bw_pic, 456, 452, 390, 29, width);

        // посмотрим на промежуточные картинки
    write_png("contrast.png", finish, width, height);
    
    // поиграли с Гауссом
    Gauss_blur(bw_pic, blr_pic, width, height); 
    // посмотрим на промежуточные картинки
    write_png("gauss.png", finish, width, height);
    
    // сделали еще что-нибудь
    // .....
    // ....
    // ....
    // ....
    // ....
    // ....
    // ....
    //
    
    write_png("intermediate_result.png", finish, width, height);
    color(blr_pic, finish, bw_size); 
    
    // выписали результат
    write_png("picture_out.png", finish, width, height); 

    printf("Всего танкеров на изображении: %d", all_tankers);
    
    // не забыли почистить память!
    free(bw_pic); 
    free(blr_pic); 
    free(finish); 
    free(picture); 
    
    return 0; 
}
