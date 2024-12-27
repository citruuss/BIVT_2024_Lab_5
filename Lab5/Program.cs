using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Security.Cryptography;
using System.Text;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public static void Main()
    {
        Program program = new Program();
    }
    #region Level 1
    public long Task_1_1(int n, int k)
    {
        long answer = 0;

        // code here
        if (n < 0 || k < 0 || n < k) return 0;
        answer = combinations(n, k);

        // create and use Combinations(n, k);
        // create and use Factorial(n);

        // end

        return answer;
    }
    public long combinations(int n, int k)
    {
        if (n < k) return 0;
        return factorial(n) / factorial(k) * factorial(n - k);
    }
    public long factorial(int n)
    {
        if (n == 0 || n == 1) return 1;
        int result = 1;
        for (int i = 2; i < n; i++)
        {
            result *= i;
        }
        return result;
    }


    public int Task_1_2(double[] first, double[] second)
    {
        int answer = 0;

        // code here
        if (first.Length != 3 || second.Length != 3) return -1;
        for (int i = 0; i < 3; i++)
        {
            if (first[i] <= 0 || second[i] <= 0)
                return -1;
        }
        double a1 = first[0], b1 = first[1], c1 = first[2];
        double a2 = second[0], b2 = second[1], c2 = second[2];
        if (a1 >= b1 + c1 || b1 >= a1 + c1 || c1 >= a1 + b1 || a2 >= b2 + c2 || b2 >= a2 + c2 || c2 >= a2 + b2) return -1;

        double Area1 = GeronArea(first[0], first[1], first[2]), Area2 = GeronArea(second[0], second[1], second[2]);
        if (Area1 == Area2) answer = 0;
        else if (Area1 > Area2) answer = 1;
        else answer = 2;
        // create and use GeronArea(a, b, c);
      
        // end

        // first = 1, second = 2, equal = 0, error = -1
        return answer;
    }
    public double GeronArea(double a, double b, double c)
    {
        double p = (a + b + c) / 2;
        double area = Math.Sqrt(p * (p - a) * (p - b) * (p - c));
        return area;
    }


    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        int answer = 0;

        // code here
        if (time <= 0 || v1 <= 0 || v2 <= 0 || a1 <= 0 || a2 <= 0)
        {
            return -1;
        }
        double s1 = GetDistance(v1, a1, time);
        double s2 = GetDistance(v2, a2, time);
        if (s1 > s2) answer = 1;
        else if (s2 > s1) answer = 2;
        else return 0;
        // create and use GetDistance(v, a, t); t - hours
        
        // end

        // first = 1, second = 2, equal = 0
        return answer;
    }
    public double GetDistance(double v, double a, double t)
    {
        double s = v * t + a * t * t / 2;
        return s;
    }

    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int answer = 0;

        // code here
        for (int time = 1; ; time++)
        {
            double s1 = GetDistance(v1, a1, time);
            double s2 = GetDistance(v2, a2, time);

            if (s1 <= s2)
            {
                return time;
            }
        }
        // use GetDistance(v, a, t); t - hours
        
        // end

        return answer;
    }
    #endregion

    #region Level 2
    public void Task_2_1(int[,] A, int[,] B)
    {
        // code here
        int amaxi, amaxj;
        int bmaxi, bmaxj;
        FindMaxIndex(A, out amaxi, out amaxj);
        FindMaxIndex(B, out bmaxi, out bmaxj);
        (A[amaxi, amaxj], B[bmaxi, bmaxj]) = (B[bmaxi, bmaxj], A[amaxi, amaxj]);

        // create and use FindMaxIndex(matrix, out row, out column);
        
        // end
    }
    public void FindMaxIndex(int[,] matrix, out int maxi, out int maxj)
    {
        int max = -100;
        maxi = -1;
        maxj = -1;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    maxi = i;
                    maxj = j;
                }
            }
        }
    }

    public void Task_2_2(double[] A, double[] B)
    {
        // code here

        // create and use FindMaxIndex(array);
        // only 1 array has to be changed!

        // end
    }

    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        // code here
        int[,] newB = new int[4, 5];
        int[,] newC = new int[5, 6];
        int Bmax = FindDiagonalMaxIndex(B);
        int Cmax = FindDiagonalMaxIndex(C);
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                if (i < Bmax)
                    newB[i, j] = B[i, j];
                else
                    newB[i, j] = B[i + 1, j];
            }
        }
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                if (i < Cmax)
                    newC[i, j] = C[i, j];
                else
                    newC[i, j] = C[i + 1, j];
            }
        }
        B = newB;
        C = newC;
        //  create and use method FindDiagonalMaxIndex(matrix);

        // end
    }
    public int FindDiagonalMaxIndex(int[,] matrix)
    {
        int max = -1000;
        int maxi = -1;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            if (max < matrix[i, i])
            {
                    max = matrix[i, i];
                    maxi = i;
            }
        }
        return maxi;
    }
    
    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here

        //  create and use method FindDiagonalMaxIndex(matrix); like in Task_2_3

        // end
    }

     public void Task_2_5(int[,] A, int[,] B)
     {
        // code here
        int Amax = FindMaxInColumn(A, 0);
        int Bmax = FindMaxInColumn(B, 0);
        for (int j = 0; j < 6; j++)
        {
            (A[Amax, j], B[Bmax, j]) = (B[Bmax, j],A[Amax, j]);
        }

        // create and use FindMaxInColumn(matrix, columnIndex, out rowIndex);
        
        // end
    }
    public int FindMaxInColumn(int[,] matrix, int columnIndex)
    {
        int max = -100;
        int maxi = -1;

        int n = matrix.GetLength(0);
        for (int i = 0; i < n; i++)
        {
            if (matrix[i, columnIndex] > max)
            {
                max = matrix[i, columnIndex];
                maxi = i;
            }
        }
        return maxi;
    }
    public void Task_2_6(ref int[] A, int[] B)
     {
        // code here

        // create and use FindMax(matrix, out row, out column); like in Task_2_1
        // create and use DeleteElement(array, index);

        // end
     }

     public void Task_2_7(ref int[,] B, int[,] C)
     {
        // code here
        int Bmax = -1;
        int Bmaxi = -1;
        for (int i = 0; i < 4; i++)
        {
            if (CountRowPositive(B, i) > Bmax)
            {
                Bmax = CountRowPositive(B, i);
                Bmaxi = i;
            }
        }
        int Cmax = -1;
        int Cmaxi = -1;
        for (int j = 0; j < 6; j++)
        {
            if (CountColumnPositive(C, j) > Cmax)
            {
                Cmax = CountColumnPositive(C, j);
                Cmaxi = j;
            }
        }
        int[,] resB = new int[5, 5];
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                if (i <= Bmaxi) resB[i, j] = B[i, j];
                else if (i == Bmaxi + 1) resB[i, j] = C[j, Cmaxi];
                else resB[i, j] = B[i - 1, j];
            }
        }
        B = resB;
        // create and use CountRowPositive(matrix, rowIndex);
        
        // create and use CountColumnPositive(matrix, colIndex);
        
        // end
    }
    public int CountRowPositive(int[,] matrix, int rowIndex)
    {
        int res = 0;
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            if (matrix[rowIndex, j] > 0)
                res++;
        }
        return res;
    }
    public int CountColumnPositive(int[,] matrix, int colIndex)
    {
        int res = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            if (matrix[i, colIndex] > 0)
                res++;
        }
        return res;
    }
    public void Task_2_8(int[] A, int[] B)
    {
        // code here

        // create and use SortArrayPart(array, startIndex);

        // end
    }

    public int[] Task_2_9(int[,] A, int[,] C)
    {
        int[] answer = default(int[]);

        // code here
        answer = new int[A.GetLength(1) + C.GetLength(1)];
        int[] s1 = SumPositiveElementsInColumns(A);
        int[] s2 = SumPositiveElementsInColumns(C);
        for (int i=0; i<s1.Length; i++)
        {
            answer[i] = s1[i];
        }
        for (int i=0; i<s2.Length; i++)
        {
            answer[i + s1.Length] = s2[i];
        }
        // create and use SumPositiveElementsInColumns(matrix);
        
        // end

        return answer;
    }
    public int[] SumPositiveElementsInColumns(int[,] matrix)
    {
        int n = matrix.GetLength(0);
        int m = matrix.GetLength(1);
        int[] ans = new int[m];
        for (int j = 0; j < m; j++)
        {
            int sum = 0;
            for (int i = 0; i < n; i++)
            {
                if (matrix[i, j] > 0)
                {
                    sum += matrix[i, j];
                }
            }
            ans[j] = sum;
        }
        return ans;
    }
    public void Task_2_10(ref int[,] matrix)
    {
        // code here

        // create and use RemoveColumn(matrix, columnIndex);

        // end
    }

    public void Task_2_11(int[,] A, int[,] B)
    {
        // code here
        FindMaxIndex(A, out int Arow, out int Acol);
        FindMaxIndex(B, out int Brow, out int Bcol);
        (A[Arow, Acol], B[Brow, Bcol]) = (B[Brow, Bcol], A[Arow, Acol]);
        // use FindMaxIndex(matrix, out row, out column); from Task_2_1
        
        // end
    }
    
    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxColumnIndex(matrix);

        // end
    }

    public void Task_2_13(ref int[,] matrix)
    {
        // code here
        int maxind = 0, minind = 0;
        int max = matrix[0, 0], min = matrix[0, 0];
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    maxind = i;
                }
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    minind = i;
                }
            }
        }
        if (maxind == minind) RemoveRow(ref matrix, maxind);
        else
        {
            if (maxind > minind)
            {
                RemoveRow(ref matrix, maxind);
                RemoveRow(ref matrix, minind);
            }
            else
            {
                RemoveRow(ref matrix, minind);
                RemoveRow(ref matrix, maxind);
            }
        }
        // create and use RemoveRow(matrix, rowIndex);
        
        // end
    }
    public void RemoveRow(ref int[,] matrix, int rowIndex)
    {
        int n = matrix.GetLength(0);
        int m = matrix.GetLength(1);
        int[,] newmatrix = new int[n - 1, m];
        int k = 0;
        for (int i = 0; i < n; i++)
        {
            if (i == rowIndex)
            {
                continue;
            }
            for (int j = 0; j < m; j++)
            {
                newmatrix[k, j] = matrix[i, j];
            }
            k++;
        }

        matrix = newmatrix;
    }

    public void Task_2_14(int[,] matrix)
    {
        // code here

        // create and use SortRow(matrix, rowIndex);

        // end
    }

    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        // code here
        double[] m = { GetAverageWithoutMinMax(A), GetAverageWithoutMinMax(B), GetAverageWithoutMinMax(C) };
        if (m[0] < m[1] && m[1] < m[2])
            answer = 1;
        else if (m[0] > m[1] && m[1] > m[2])
            answer = -1;
        else
            answer = 0;
        // create and use GetAverageWithoutMinMax(matrix);
       
        // end

        // 1 - increasing   0 - no sequence   -1 - decreasing
        return answer;
    }
     public double GetAverageWithoutMinMax(int[,] matrix)
    {
        int n = matrix.GetLength(0);
        int m = matrix.GetLength(1);
        int total = n * m;
        if (total <= 2)
        {
            return 0;
        }
        int min = 1000, max = -1000, sum = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                int p = matrix[i, j];
                sum += p;
                if (p < min) min = p;
                if (p > max) max = p;
            }
        }
        sum -= (min + max);
        int answer = sum / (total - 2);
        return answer;
    }
    public void Task_2_16(int[] A, int[] B)
    {
        // code here

        // create and use SortNegative(array);

        // end
    }

    public void Task_2_17(int[,] A, int[,] B)
    {
        // code here
        SortRowsByMaxElement(A);
        SortRowsByMaxElement(B);
        // create and use SortRowsByMaxElement(matrix);

        // end
    }
    public void SortRowsByMaxElement(int[,] matrix)
    {
        for (int i = 1; i < matrix.GetLength(0); i++)
        {
            for (int j = i; j > 0; j--)
            {
                if (matrix[j - 1, FindRowMaxIndex(matrix, j - 1)] < matrix[j, FindRowMaxIndex(matrix, j)])
                {
                    for (int col = 0; col < matrix.GetLength(1); col++)
                    {
                        (matrix[j - 1, col], matrix[j, col]) = (matrix[j, col], matrix[j - 1, col]);
                    }
                }
                else
                    break;
            }
        }
    }
    public int FindRowMaxIndex(int[,] matrix, int rowIndex)
    {
        int max = -1000;
        int maxi = -1;
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            if (matrix[rowIndex, j] > max)
            {
                max = matrix[rowIndex, j];
                maxi = j;
            }
        }
        return maxi;
    }
    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here

        // create and use SortDiagonal(matrix);

        // end
    }

    public void Task_2_19(ref int[,] matrix)
    {
        // code here
        if (matrix == null || matrix.GetLength(0) == 0) 
        { 
            return; 
        }
        for (int i = matrix.GetLength(0) - 1; i >= 0; i--)
        {
            int k = 0;
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] == 0)
                {
                    k = 1;
                    break;
                }
            }
            if (k == 1)
            {
                RemoveRow(ref matrix, i);
            }
        }

        // use RemoveRow(matrix, rowIndex); from 2_13

        // end
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        // use RemoveColumn(matrix, columnIndex); from 2_10

        // end
    }

    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = null;
        answerB = null;

        // code here
        answerA = CreateArrayFromMins(A);
        answerB = CreateArrayFromMins(B);
        
        // create and use CreateArrayFromMins(matrix);

        // end
    }
    public int[] CreateArrayFromMins(int[,] matrix)
    {
        int[] res = new int[matrix.GetLength(0)];
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int min = 1000;
            for (int j = i; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                }
                res[i] = min;
            }
        }
        return res;
    }
    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here

        // create and use CountNegativeInRow(matrix, rowIndex);
        // create and use FindMaxNegativePerColumn(matrix);

        // end
    }

    public void Task_2_23(double[,] A, double[,] B)
    {
        // code here
        MatrixValuesChange(A);
        MatrixValuesChange(B);
        // create and use MatrixValuesChange(matrix);

        // end
    }
    public void MatrixValuesChange(double[,] matrix)
    {
        int n = matrix.GetLength(0);
        int m = matrix.GetLength(1);
        double[] array = new double[n * m];
        int index = 0;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                array[index++] = matrix[i, j];
            }
        }
        for (int i = 0; i < array.Length; i++)
        {
            for (int j = 0; j < array.Length - i - 1; j++)
            {
                if (array[j] < array[j + 1])
                {
                    double temp = array[j];
                    array[j] = array[j + 1];
                    array[j + 1] = temp;
                }
            }
        }
        if (array.Length < 5)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    if (matrix[i, j] > 0) 
                    { 
                        matrix[i, j] *= 2; 
                    }
                    else matrix[i, j] *= 0.5;
                }
            }
        }
        else
        {
            double[] max = new double[5];
            for (int i = 0; i < 5; i++)
            {
                max[i] = array[i];
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    bool flag = false;
                    for (int k = 0; k < max.Length; k++)
                    {
                        if (matrix[i, j] == max[k])
                        {
                            flag = true;
                            break;
                        }
                    }
                    if (flag)
                    {
                        if (matrix[i, j] > 0) 
                        { 
                            matrix[i, j] *= 2; 
                        }
                        else matrix[i, j] *= 0.5;
                    }
                    else
                    {
                        if (matrix[i, j] < 0) 
                        { 
                            matrix[i, j] *= 2; 
                        }
                        else matrix[i, j] /= 2;
                    }
                }
            }
        }
    }
    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); like in 2_1
        // create and use SwapColumnDiagonal(matrix, columnIndex);

        // end
    }

    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = 0;
        maxB = 0;

        // code here
        maxA = FindRowWithMaxNegativeCount(A);
        maxB = FindRowWithMaxNegativeCount(B);
        // create and use FindRowWithMaxNegativeCount(matrix);
        // in FindRowWithMaxNegativeCount create and use CountNegativeInRow(matrix, rowIndex); like in 2_22

        // end
    }
    public int FindRowWithMaxNegativeCount(int[,] matrix)
    {
        int max = -1000;
        int row = 0;
        for(int i=0; i<matrix.GetLength(0); i++) 
        { 
            if (CountNegativeInRow(matrix, i) > max)
            {
                max = CountNegativeInRow(matrix, i);
                row = i;
            }
        }
        return row;
        
    }
    public int CountNegativeInRow(int[,] matrix, int rowIndex)
    {
        int k = 0;
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            if (matrix [rowIndex, j] < 0)
            {
                k++;
            }
        }
        return k;
    }
    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22

        // end
    }

    public void Task_2_27(int[,] A, int[,] B)
    {
        // code here
        ReplaceMaxElementEven(A);
        ReplaceMaxElementOdd(A);
        ReplaceMaxElementEven(B);
        ReplaceMaxElementOdd(B);
        // create and use FindRowMaxIndex(matrix, rowIndex, out columnIndex);
        // create and use ReplaceMaxElementOdd(matrix, row, column);
        // create and use ReplaceMaxElementEven(matrix, row, column);

        // end
    }
    public void ReplaceMaxElementEven(int[,] matrix)
    {
        int row = -1;
        for (int i = 1; i < matrix.GetLength(0); i += 2)
        {
            row = FindRowMaxIndex(matrix, i);
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] == matrix[i, row])
                {
                    matrix[i, j] = 0;
                }
            }
        }
    }
    public void ReplaceMaxElementOdd(int[,] matrix)
    {
        int row = -1;
        for (int i = 0; i < matrix.GetLength(0); i += 2)
        {
            row = FindRowMaxIndex(matrix, i);
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] == matrix[i, row])
                {
                    matrix[i, j] *= j + 1;
                }
                    
            }
        }
    }
    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create and use FindSequence(array, A, B); // 1 - increasing, 0 - no sequence,  -1 - decreasing
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a or Task_2_28b
        // A and B - start and end indexes of elements from array for search

        // end
    }
    #endregion

    #region Level 3
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        // code here
        double a1 = 0.1, b1 = 1, h1 = 0.1;
        firstSumAndY = new double[(int)((b1 - a1) / h1) + 1, 2];
        GetSumAndY(S1, y1, a1, b1, h1, firstSumAndY);
        double a2 = Math.PI / 5, b2 = Math.PI, h2 = Math.PI / 25;
        secondSumAndY = new double[(int)((b2 - a2) / h2) + 1, 2];
        GetSumAndY(S2, y2, a2, b2, h2, secondSumAndY, 1);

        // create and use public delegate SumFunction(x) and public delegate YFunction(x)
        // create and use method GetSumAndY(sFunction, yFunction, a, b, h);
        // create and use 2 methods for both functions calculating at specific x

        // end
    }
    public delegate double SumFunction(int i, double x, ref int change);
    public delegate double YFunction(double x);
    public double S1(int i, double x, ref int iFact)
    {
        if (i > 0)
        {
            iFact *= i;
        }
        return Math.Cos(i * x) / iFact;
    }
    public double S2(int i, double x, ref int sign)
    {
        sign *= -1;
        return sign * Math.Cos(i * x) / (i * i);
    }
    public double y1(double x)
    {
        return Math.Exp(Math.Cos(x)) * Math.Cos(Math.Sin(x));
    }
    public double y2(double x)
    {
        return ((x * x) - Math.PI * Math.PI / 3) / 4;
    }
    public double CalculateSum(SumFunction sumFunction, double x, int i)
    {
        double eps = 0.0001, sum = 0;
        int ch = 1;
        double c = sumFunction(i, x, ref ch);
        while (Math.Abs(c) > eps)
        {
            sum += c;
            c = sumFunction(++i, x, ref ch);
        }
        return sum;
    }
    public void GetSumAndY(SumFunction sFunction, YFunction yFunction, double a, double b, double h, double[,] SumAndY, int startI = 0)
    {
        for (int i = 0; i < (b - a) / h + 1; i++)
        {
            double x = a + i * h;
            double sum = CalculateSum(sFunction, x, startI);
            double y = yFunction(x);
            SumAndY[i, 0] = sum;
            SumAndY[i, 1] = y;
        }
    }

    public void Task_3_2(int[,] matrix)
    {
        // SortRowStyle sortStyle = default(SortRowStyle); - uncomment me

        // code here

        // create and use public delegate SortRowStyle(matrix, rowIndex);
        // create and use methods SortAscending(matrix, rowIndex) and SortDescending(matrix, rowIndex)
        // change method in variable sortStyle in the loop here and use it for row sorting

        // end
    }

    public double Task_3_3(double[] array)
    {
        double answer = 0;
        SwapDirection swapper = default(SwapDirection);
        // code here
        if (array == null || array.Length == 0) 
        { 
            return 0; 
        }
        double average = CalculateAverage(array);
        if (array[0] > average) 
        { 
            swapper = SwapLeft; 
        }
        else swapper = SwapRight;
        swapper(array);
        answer = GetSum(array);
        // create and use public delegate SwapDirection(array);
        // create and use methods SwapRight(array) and SwapLeft(array)
        // create and use method GetSum(array, start, h) that sum elements with a negative indexes
        // change method in variable swapper in the if/else and than use swapper(matrix)

        // end

        return answer;
    }
    public delegate void SwapDirection(double[] array);
    public double CalculateAverage(double[] array)
    {
        double sum = 0;
        foreach (int n in array)
        {
            sum += n;
        }
        return sum / array.Length;
    }
    public void SwapLeft(double[] array)
    {
        for (int i = 0; i < array.Length - 1; i += 2)
        {
            double temp = array[i];
            array[i] = array[i + 1];
            array[i + 1] = temp;
        }
    }
    public void SwapRight(double[] array)
    {
        for (int i = array.Length - 1; i > 0; i -= 2)
        {
            double temp = array[i];
            array[i] = array[i - 1];
            array[i - 1] = temp;
        }
    }
    public double GetSum(double[] array)
    {
        double sum = 0;
        for (int i = 1; i < array.Length; i += 2)
        {
            sum += array[i];
        }
        return sum;
    }

    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here

        // create and use public delegate GetTriangle(matrix);
        // create and use methods GetUpperTriange(array) and GetLowerTriange(array)
        // create and use GetSum(GetTriangle, matrix)

        // end

        return answer;
    }

    public void Task_3_5(out int func1, out int func2)
    {
        func1 = 0;
        func2 = 0;

        // code here
        double a1 = 0, b1 = 2, h1 = 0.1;
        double a2 = -1, b2 = 1, h2 = 0.2;
        func1 = CountSignFlips(F1, a1, b1, h1);
        func2 = CountSignFlips(F2, a2, b2, h2);
        // use public delegate YFunction(x, a, b, h) from Task_3_1
        // create and use method CountSignFlips(YFunction, a, b, h);
        // create and use 2 methods for both functions

        // end
    }
    public delegate double FunctionDelegate(double x);
    public double F1(double x)
    {
        return x * x - Math.Sin(x);
    }
    public double F2(double x)
    {
        return Math.Exp(x) - 1;
    }
    public int CountSignFlips(YFunction y, double a, double b, double h)
    {
        int count = 0;
        double last = y(a);
        for (double x = a + h; x <= b; x += h)
        {
            double pr = y(x);
            if (last * pr < 0)
            {
                count++;
            }
            if (pr != 0)
            {
                last = pr;
            }
        }
        return count;
    }
    public void Task_3_6(int[,] matrix)
    {
        // code here

        // create and use public delegate FindElementDelegate(matrix);
        // use method FindDiagonalMaxIndex(matrix) like in Task_2_3;
        // create and use method FindFirstRowMaxIndex(matrix);
        // create and use method SwapColumns(matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);

        // end
    }

    public delegate int CountPositive(int[,] matrix, int index);
    public void Task_3_7(ref int[,] B, int[,] C)
    {
        // code here
        CountPositive count = CountRowPositive;
        int Bmax = -1;
        int Bmaxi = -1;
        for (int i = 0; i < 4; i++)
        {
            if (count(B, i) > Bmax)
            {
                Bmax = count(B, i);
                Bmaxi = i;
            }
        }
        count = CountColumnPositive;
        int Cmax = -1;
        int Cmaxi = -1;
        for (int j = 0; j < 6; j++)
        {
            if (count(C, j) > Cmax)
            {
                Cmax = count(C, j);
                Cmaxi = j;
            }
        }
        int[,] Bres = new int[5, 5];
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                if (i <= Bmaxi)
                    Bres[i, j] = B[i, j];
                else if (i == Bmaxi + 1)
                    Bres[i, j] = C[j, Cmaxi];
                else
                    Bres[i, j] = B[i - 1, j];
            }
        }
        B = Bres;
        // create and use public delegate CountPositive(matrix, index);
        // use CountRowPositive(matrix, rowIndex) from Task_2_7
        // use CountColumnPositive(matrix, colIndex) from Task_2_7
        // create and use method InsertColumn(matrixB, CountRow, matrixC, CountColumn);

        // end
    }

    public void Task_3_10(ref int[,] matrix)
    {
        // FindIndex searchArea = default(FindIndex); - uncomment me

        // code here

        // create and use public delegate FindIndex(matrix);
        // create and use method FindMaxBelowDiagonalIndex(matrix);
        // create and use method FindMinAboveDiagonalIndex(matrix);
        // use RemoveColumn(matrix, columnIndex) from Task_2_10
        // create and use method RemoveColumns(matrix, findMaxBelowDiagonalIndex, findMinAboveDiagonalIndex)

        // end
    }

    public void Task_3_13(ref int[,] matrix)
    {
        // code here
        RemoveRows(ref matrix, FindMax, FindMin);
        // use public delegate FindElementDelegate(matrix) from Task_3_6
        // create and use method RemoveRows(matrix, findMax, findMin)

        // end
    }
    public delegate int FindElementDelegate(int[,] matrix);
    public int FindMax(int[,] matrix)
    {
        int max = matrix[0, 0];
        int maxi = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    maxi = i;
                }
            }
        }
        return maxi;
    }
    public int FindMin(int[,] matrix)
    {
        int min = matrix[0, 0];
        int mini = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    mini = i;
                }
            }
        }
        return mini;
    }
    public void RemoveRows(ref int[,] matrix, FindElementDelegate FindMax, FindElementDelegate FindMin)
    {
        int maxi = FindMax(matrix);
        int mini = FindMin(matrix);
        if (maxi == mini)
        {
            RemoveRow(ref matrix, maxi);
        }
        else
        {
            RemoveRow(ref matrix, Math.Max(maxi, mini));
            RemoveRow(ref matrix, Math.Min(maxi, mini));
        }
    }
    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {

        rows = null;
        cols = null;

        // code here

        // create and use public delegate GetNegativeArray(matrix);
        // use GetNegativeCountPerRow(matrix) from Task_2_22
        // use GetMaxNegativePerColumn(matrix) from Task_2_22
        // create and use method FindNegatives(matrix, searcherRows, searcherCols, out rows, out cols);

        // end
    }

    public void Task_3_27(int[,] A, int[,] B)
    {
        // code here
        EvenOddRowsTransform(A, ReplaceMaxElementOdd, ReplaceMaxElementEven);
        EvenOddRowsTransform(B, ReplaceMaxElementOdd, ReplaceMaxElementEven);
        // create and use public delegate ReplaceMaxElement(matrix, rowIndex, max);
        // use ReplaceMaxElementOdd(matrix) from Task_2_27
        // use ReplaceMaxElementEven(matrix) from Task_2_27
        // create and use method EvenOddRowsTransform(matrix, replaceMaxElementOdd, replaceMaxElementEven);

        // end
    }
    public delegate void ReplaceMaxElement(int[,] matrix);
    public void EvenOddRowsTransform(int[,] matrix, ReplaceMaxElement replaceMaxElementEven, ReplaceMaxElement replaceMaxElementOdd)
    {
        replaceMaxElementEven(matrix);
        replaceMaxElementOdd(matrix);
    }
    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // create and use method FindIncreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method FindDecreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method DefineSequence(array, findIncreasingSequence, findDecreasingSequence);

        // end
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // use method FindIncreasingSequence(array, A, B); from Task_3_28a
        // use method FindDecreasingSequence(array, A, B); from Task_3_28a
        // create and use method FindLongestSequence(array, sequence);

        // end
    }
    #endregion
    #region bonus part
    public double[,] Task_4(double[,] matrix, int index)
    {
        if (matrix.GetLength(0) != matrix.GetLength(1)) return matrix;
        MatrixConverter[] mc = new MatrixConverter[4] { ToUpperTriangular, ToLowerTriangular, ToLeftDiagonal, ToRightDiagonal };
        mc[index](matrix);
        return matrix;
        // code here

        // create public delegate MatrixConverter(matrix);
        // create and use method ToUpperTriangular(matrix);
        // create and use method ToLowerTriangular(matrix);
        // create and use method ToLeftDiagonal(matrix); - start from the left top angle
        // create and use method ToRightDiagonal(matrix); - start from the right bottom angle

        // end
    }
    public delegate void MatrixConverter(double[,] matrix);
    public void ToUpperTriangular(double[,] matrix)
    {
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            for (int i = j + 1; i < matrix.GetLength(0); i++)
            {
                double p = -(matrix[i, j] / matrix[j, j]);
                matrix[i, j] = 0;
                for (int k = j + 1; k < matrix.GetLength(1); k++)
                {
                    matrix[i, k] += matrix[j, k] * p;
                }
            }
        }
    }
    public void ToLowerTriangular(double[,] matrix)
    {
        for (int j = matrix.GetLength(1) - 1; j >= 0; j--)
        {
            for (int i = j - 1; i >= 0; i--)
            {
                double p = -(matrix[i, j] / matrix[j, j]);
                matrix[i, j] = 0;
                for (int k = j - 1; k >= 0; k--)
                {
                    matrix[i, k] += matrix[j, k] * p;
                }
            }
        }
    }
    public void ToLeftDiagonal(double[,] matrix)
    {
        ToUpperTriangular(matrix);
        ToLowerTriangular(matrix);
    }
    public void ToRightDiagonal(double[,] matrix)
    {
        ToLowerTriangular(matrix);
        ToUpperTriangular(matrix);
    }
    #endregion
}
