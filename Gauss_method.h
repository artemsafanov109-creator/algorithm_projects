#include <stdio.h>

void addition(int len, double *matrix, int index1, int index2, int flag) {
    for (int i = 0; i<len; i++){
        if (flag) {matrix[i*len + index1] = matrix[i*len + index1] + matrix[i*len + index2];}
        else {matrix[index1*len + i] = matrix[index1*len + i] + matrix[index2*len + i];}
    }
}

void multiplication(double value, int len, double *matrix, int index, int flag) {
    for (int i = 0; i<len; i++){
        if (flag) {matrix[i*len + index] = matrix[i*len + index] * value;}
        else {matrix[index*len + i] = matrix[index*len + i] * value;}
    }
}

void change_lines(int len, double *matrix, int index1, int index2, int flag) {
    for (int i = 0; i<len; i++){
        int t;
        if (flag) {
            t = matrix[i*len + index2];
            matrix[i*len + index2] = matrix[i*len + index1];
            matrix[i*len + index1] = t;
            
        }
        else {
            t = matrix[index2*len + i];
            matrix[index2*len + i] = matrix[index1*len + i];
            matrix[index1*len + i] = t;
        }
    }
}

void Gauss_Method(int len, double *matrix, double *det_change) {
    if (len == 1) {return;}
    int i = 0, j = 0;
    double det_main = *det_change;
    while (matrix[i*len + j] == 0) {
        i++;
        if (i == len) {j++; i = 0;}
        if (j == len) {return;}
    }
    if (i != 0) {change_lines(len, matrix, i, 0, 0); det_main = det_main * (-1);}
    while (i != len - 1) {
        i++;
        if (matrix[i*len + j] != 0) {
            double lambda = (-1)*(matrix[j] / matrix[i*len + j]);
            multiplication(lambda, len, matrix, i, 0);
            det_main = det_main * lambda;
            addition(len, matrix, i, 0, 0);
        }
    }
    
    int size = (len - j - 1) * (len - j - 1); // Создадим доп матрицу, в которую перепишем матрицу [a_ij], где i = 2,...,len , а j = j+1,..., len-1 (в нашем случае i = 1,...,len-1, j = j+1,..., len-1)
    double extra_matrix[size];
    int j1 = j;
    int k = 0; // индекс в новой матрице
    for (int w=1; w<len; w++) {
        for (int t=j1+1; t<len; t++) {
            extra_matrix[k++] = matrix[w*len + t];
        }
    }
    Gauss_Method(len - 1, extra_matrix, &det_main);
    k = 0; // тот же индекс k
    for (int w=1; w<len; w++) {
        for (int t=j1+1; t<len; t++) {
            matrix[w*len + t] = extra_matrix[k++];
        }
    }
    *det_change = det_main;
    return;
}

int main() {
    int hight;
    scanf("%d", &hight);
    int matrsize = hight * hight;
    double main_matrix[matrsize];
    for (int i=0; i<matrsize; i++) {
        scanf("%lf", &main_matrix[i]);
    }
    double det_mark = 1;
    Gauss_Method(hight, main_matrix, &det_mark);
    double det = 1.0;
    for (int i=0; i<hight; i++) {
        if (main_matrix[i*hight + i] == 0.0 || main_matrix[i*hight + i] == -0.0) {det = 0.0;}
        det *= main_matrix[i*hight + i];
    }
    det = det / det_mark;
    printf("%lf", det);
    return 0;
}
