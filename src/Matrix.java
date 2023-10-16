import java.io.*;
import java.util.*;
import java.text.*;
public class Matrix {
    protected int rowsCount;
    protected int columnsCount;
    protected double[][] matrix;
    Matrix(String pathToFile) throws FileNotFoundException
    {
        File input = new File(pathToFile);
        Scanner whileScan = new Scanner(input);
        String line = whileScan.nextLine();
        String[] strArr = line.trim().split("\\s+");
        double[] row;

        this.rowsCount = 1;
        this.columnsCount = strArr.length;
        this.matrix = new double[this.rowsCount][this.columnsCount];


        while (whileScan.hasNextLine())
        {
            line = whileScan.nextLine();
            strArr = line.trim().split("\\s+");
            row = new double[strArr.length];
            for (int i = 0; i < strArr.length; i++)
                row[i] = Double.parseDouble(strArr[i]);
            this.addRow(row);
        }
        Scanner finalScan = new Scanner(input);
        line = finalScan.nextLine();
        strArr = line.trim().split("\\s+");
        row = new double[strArr.length];
        for (int i = 0; i < strArr.length; i++)
            row[i] = Double.parseDouble(strArr[i]);
        this.setRow(row, 0);
    }
    Matrix()
    {
        System.out.println(Main.INPUT + "Введите количество строк и столбцов в матрице:" + Main.RESET);
        Scanner scan = new Scanner(System.in);
        int rowsCount = scan.nextInt();
        int columnsCount = scan.nextInt();
        this.rowsCount = rowsCount;
        this.columnsCount = columnsCount;
        this.matrix = new double[rowsCount][columnsCount];
    }
    Matrix(int rowsCount, int columnsCount)
    {
        this.rowsCount = rowsCount;
        this.columnsCount = columnsCount;
        this.matrix = new double[rowsCount][columnsCount];
    }
    Matrix(double[][] matrix, int rowsCount, int columnsCount)
    {
        this.rowsCount = rowsCount;
        this.columnsCount = columnsCount;
        this.matrix = matrix;
    }
    int getRowsCount()
    { return this.rowsCount; }
    int getColumnsCount()
    { return this.columnsCount; }
    double[][] getMatrix()
    { return this.matrix; }
    void setMatrix(double[][] matrix)
    { this.matrix = matrix; }
    void setItem(double replaceItem, int rowIndex, int colIndex)
    { this.matrix[rowIndex][colIndex] = replaceItem; }
    void setRow(double[] row, int index)
    {
        if (this.columnsCount != row.length)
        { System.out.println(Main.ERROR + "Невозможно строку строку в исходной матрице из-за несоответсвия количества столбцов" + Main.RESET); }
        else
            this.matrix[index] = row;
    }
    void addRow(double[] row)
    {
        if (row.length != this.columnsCount)
        { System.out.println(Main.ERROR + "Невозможно добавить строку в исходную матрицу из-за несоответсвия количества столбцов" + Main.RESET); }
        else
        {
            this.rowsCount ++;
            double[][] newMatrix = new double[this.rowsCount][this.columnsCount];
            for (int i = 0; i < this.rowsCount; i++)
            {
                if (i != this.rowsCount - 1)
                    System.arraycopy(this.matrix[i], 0, newMatrix[i], 0, this.columnsCount);
                else newMatrix[i] = row;
            }
            this.matrix = newMatrix;
        }
    }
    void deleteRow(int index)
    {
        if (index > this.rowsCount || index < 0)
        { System.out.println(Main.ERROR + "Номер строки указан неверно" + Main.RESET); }
        else
        {
            this.rowsCount --;
            double[][] newMatrix = new double[this.rowsCount][this.columnsCount];
            for (int i = 0; i < this.rowsCount + 1; i++)
                if (i != index)


                    this.matrix = newMatrix;
        }
    }
    void addColumn(double[] column)
    {
        if (column.length != this.rowsCount)
        { System.out.println(Main.ERROR + "Невозможно добавить столбец в исходную матрицу из-за несоответсвия количества строк" + Main.RESET); }
        else
        {
            this.columnsCount ++;
            double[][] newMatrix = new double[this.rowsCount][this.columnsCount];
            for (int i = 0; i < this.rowsCount; i++)
                for (int j = 0; j < this.columnsCount; j++)
                {
                    if (j != this.columnsCount - 1)
                        newMatrix[i][j] = this.matrix[i][j];
                    else newMatrix[i][j] = column[i];
                }
            this.matrix = newMatrix;
        }
    }
    void createMatrix()
    {
        Scanner scan = new Scanner(System.in);
        System.out.println(Main.INPUT + "Введите элементы матрицы " + rowsCount + " на " + columnsCount + ":" + Main.RESET);
        for(int i = 0; i < this.rowsCount; i++)
            for (int j = 0; j < this.columnsCount; j++)
            {
                double a = scan.nextDouble();
                this.matrix[i][j] = a;
            }
    }
    void createRandomMatrix(double from, double to)
    {
        Random random = new Random();
        for (int i = 0; i < this.rowsCount; i++)
            for (int j = 0; j < this.columnsCount; j++)
                this.matrix[i][j] = random.nextDouble(from, to);
    }
    void createRandomIntMatrix(int from, int to)
    {
        Random random = new Random();
        for (int i = 0; i < this.rowsCount; i++)
            for (int j = 0; j < this.columnsCount; j++)
                this.matrix[i][j] = random.nextInt(from, to);
    }
    void printMatrix()
    {
        System.out.println(Main.HEADER_OUTPUT + "Матрица " + rowsCount + " на " + columnsCount + ":" + Main.OUTPUT);
        for(int i = 0; i < this.rowsCount; i++)
        {
            for (int j = 0; j < this.columnsCount; j++)
            {
                System.out.print(this.matrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.print(Main.RESET);
    }
    void printFormattedMatrix()
    {
        System.out.println(Main.HEADER_OUTPUT + "Матрица " + rowsCount + " на " + columnsCount + ":" + Main.OUTPUT);
        for(int i = 0; i < this.rowsCount; i++)
        {
            for (int j = 0; j < this.columnsCount; j++)
            {
                DecimalFormat shortOut = new DecimalFormat("#.##");
                String result = shortOut.format(this.matrix[i][j]);
                System.out.print(result + " ");
            }
            System.out.println();
        }
        System.out.print(Main.RESET);
    }
    Matrix matrixAddition(Matrix addMatrix)
    {
        if (this.rowsCount != addMatrix.rowsCount || this.columnsCount != addMatrix.columnsCount)
        {
            System.out.println(Main.ERROR + "Размеры матриц разные \n" + Main.COMMENT + "Пожалуйста, введите матрицы одного размера" + Main.RESET);
            return null;
        }
        else
        {
            Matrix resultMatrix = new Matrix(this.rowsCount, addMatrix.columnsCount);
            for(int i = 0; i < resultMatrix.getRowsCount(); i++)
                for(int j = 0; j < resultMatrix.getColumnsCount(); j++)
                { resultMatrix.matrix[i][j] = this.matrix[i][j] + addMatrix.matrix[i][j]; }
            return resultMatrix;
        }
    }
    Matrix constantMultiplication(double constant)
    {
        double[][] newMatrix = this.matrix;
        for (int i = 0; i < this.rowsCount; i++)
            for (int j = 0; j < this.columnsCount; j++)
            { newMatrix[i][j] *= constant; }
        return new Matrix(newMatrix, this.rowsCount, this.columnsCount);
    }
    Matrix matrixMultiplication(Matrix addMatrix)
    {
        if (this.columnsCount != addMatrix.getRowsCount())
        {
            System.out.println(Main.ERROR + "Размеры матриц разные \n" + Main.COMMENT +
                    "Пожалуйста, убедитесь, что количество столбцов первой матрицы равно количеству строк второй матрицы" + Main.RESET);
            return null;
        }
        else
        {
            Matrix resultMatrix = new Matrix(this.rowsCount, addMatrix.getColumnsCount());
            for (int i = 0; i < this.rowsCount; i++)
            {
                resultMatrix.matrix[i] = new double[addMatrix.getColumnsCount()];
                for (int j = 0; j < addMatrix.getColumnsCount(); j++)
                {
                    resultMatrix.matrix[i][j] = 0;

                    for (int k = 0; k < addMatrix.getRowsCount(); k++)
                    { resultMatrix.matrix[i][j] += this.matrix[i][k] * addMatrix.matrix[k][j]; }
                }
            }
            return resultMatrix;
        }
    }
    Vector matrixAndVectorMultiplication(Vector addVector)
    {
        Matrix vectorMatrix = addVector.vectorToMatrix();
        Vector resultVector;
        resultVector = this.matrixMultiplication(vectorMatrix).matrixToVector();
        return resultVector;
    }
    Vector matrixToVector()
    {
        if (this.columnsCount != 1)
        {
            System.out.println(Main.ERROR + "Преобразование из матрицы в вектор невозможно." + Main.RESET);
            return null;
        }
        else
        {
            double[] convertVector = new double[this.rowsCount];
            for (int i = 0; i < this.rowsCount; i++)
                convertVector[i] = this.matrix[i][0];
            return new Vector(convertVector, this.rowsCount);
        }
    }
    Matrix transposition()
    {
        double[][] transpositionMatrix = new double[this.columnsCount][this.rowsCount];
        for(int i = 0; i < this.columnsCount; i++)
            for (int j = 0; j < this.columnsCount; j++)
            { transpositionMatrix[i][j] = this.matrix[j][i]; }
        return new Matrix(transpositionMatrix, this.columnsCount, this.rowsCount);
    }

    boolean isMatrixEqual(Matrix compareMatrix)
    {
        if (this.rowsCount != compareMatrix.getRowsCount() || this.columnsCount != compareMatrix.getColumnsCount())
        {
            System.out.println(Main.ERROR + "Размеры матриц разные" + Main.RESET);
            return false;
        }
        else
        {
            for (int i = 0; i < this.rowsCount; i++)
                for (int j = 0; j < this.columnsCount; j++)
                    if (this.matrix[i][j] != compareMatrix.getMatrix()[i][j])
                        return false;
            return true;
        }
    }
    void swapRow(int indexChange, int indexChangeWith)
    {
        double[] tmpRow = this.matrix[indexChange];
        this.matrix[indexChange] = this.matrix[indexChangeWith];
        this.matrix[indexChangeWith] = tmpRow;
    }
    void swapColumn(int indexChange, int indexChangeWith)
    {
        for (int i = 0; i < this.rowsCount; i++)
        {
            double tmpCol = this.matrix[i][indexChange];
            this.matrix[i][indexChange] = this.matrix[i][indexChangeWith];
            this.matrix[i][indexChangeWith] = tmpCol;
        }
    }
    double gaussianDeterminant()
    {
        double[][] copyMatrix = new double[this.rowsCount][this.columnsCount];
        for (int i = 0; i < this.rowsCount; i++)
            System.arraycopy(this.matrix[i], 0, copyMatrix[i], 0, this.columnsCount);
        Matrix tempMatrix = new Matrix(copyMatrix, this.rowsCount, this.columnsCount);
        double determinant = 1;
        for (int k = 0; k < rowsCount - 1; k++)
            for (int i = k + 1; i < rowsCount; i++)
            {
                if (tempMatrix.matrix[k][k] == 0)
                {
                    for (int l = i; l < this.rowsCount; l++)
                    {
                        if (tempMatrix.matrix[l][l] != 0) {tempMatrix.swapRow(k, l); determinant *= -1;}
                        else System.out.println(Main.ERROR + "Ошибка! Деление на ноль, невозможно посчитать определитель" + Main.RESET);
                    }
                }
                double tmp = tempMatrix.matrix[i][k] / tempMatrix.matrix[k][k];
                for (int j = k; j < columnsCount; j++)
                { tempMatrix.matrix[i][j] -= tmp * tempMatrix.matrix[k][j]; }
            }
        for (int i = 0; i < this.rowsCount; i++)
        { determinant *= tempMatrix.matrix[i][i]; }
        return Math.round(determinant);
    }
}