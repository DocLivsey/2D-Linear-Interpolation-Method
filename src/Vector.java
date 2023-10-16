import java.io.*;
import java.text.*;
import java.util.*;
public class Vector {
    protected double[] vector;
    protected int vectorSize;
    Vector(String pathToFile) throws FileNotFoundException
    {
        File input = new File(pathToFile);
        Scanner scan = new Scanner(input);
        String line = scan.nextLine();
        String[] strArr = line.trim().split("\\s+");

        this.vectorSize = strArr.length;
        this.vector = new double[this.vectorSize];
        for (int i = 0; i < strArr.length; i++)
            this.vector[i] = Double.parseDouble(strArr[i]);
    }
    Vector(double[] vector, int vectorSize)
    {
        this.vector = vector;
        this.vectorSize = vectorSize;
    }
    Vector(int vectorSize)
    {
        this.vectorSize = vectorSize;
        this.vector = new double[vectorSize];
    }
    Vector()
    {
        Scanner scan = new Scanner(System.in);
        System.out.println(Main.INPUT + "Введите размер вектора:" + Main.RESET);
        int vectorSize = scan.nextInt();

        this.vectorSize = vectorSize;
        this.vector = new double[vectorSize];
    }
    void createVector()
    {
        Scanner scan = new Scanner(System.in);
        System.out.println(Main.INPUT + "Введите элементы вектора размерностью " + vectorSize + ":" + Main.RESET);
        for (int i = 0; i < vectorSize; i++)
        {
            int a = scan.nextInt();
            this.vector[i] = a;
        }
    }
    void createRandomVector(double from, double to)
    {
        Random random = new Random();
        for (int i = 0; i < this.vectorSize; i++)
            this.vector[i] = random.nextDouble(from, to);
    }
    void createRandomIntVector(int from, int to)
    {
        Random random = new Random();
        for (int i = 0; i < this.vectorSize; i++)
            this.vector[i] = random.nextInt(from, to);
    }
    int getVectorSize() { return this.vectorSize; }
    double[] getVector() { return this.vector; }
    void setItem(double replaceItem, int index)
    { this.vector[index] = replaceItem; }
    void addItem(double item)
    {
        this.vectorSize ++;
        double[] newVector = new double[this.vectorSize];
        for (int i = 0; i < this.vectorSize; i++)
        {
            if (i != this.vectorSize - 1)
                newVector[i] = this.vector[i];
            else
                newVector[i] = item;
        }
        this.vector = newVector;
    }
    void setVector(double[] vector)
    {
        System.out.println(Main.CHOOSE + "Вы уверены, что хотите заменить вектор?" + Main.RESET);
        this.vector = vector;
    }
    void printVector()
    {
        System.out.print(Main.HEADER_OUTPUT + "\nВектор размерностью " + vectorSize + Main.OUTPUT + ": \n { ");
        for (int i = 0; i < vectorSize; i++)
        {
            System.out.print(this.vector[i] + "; ");
        }
        System.out.println("}" + Main.RESET);
    }
    void printFormattedVector()
    {
        System.out.print(Main.HEADER_OUTPUT + "\nВектор размерностью " + vectorSize + Main.OUTPUT + ": \n { ");
        for (int i = 0; i < vectorSize; i++)
        {
            DecimalFormat shortOut = new DecimalFormat("#.##");
            String result = shortOut.format(this.vector[i]);
            System.out.print(result + "; ");
        }
        System.out.println("}" + Main.RESET);
    }
    Vector constantMultiplication(double constant)
    {
        double[] newVector = this.vector;
        for (int i = 0; i < this.vectorSize; i++)
            newVector[i] *= constant;
        return new Vector(newVector,this.vectorSize);
    }
    Vector vectorAddition(Vector addVector)
    {
        if (this.vectorSize != addVector.getVectorSize())
        {
            System.out.println(Main.ERROR + "Размеры векторов разные \n" + Main.COMMENT + "Пожалуйста, введите вектора одного размера" + Main.RESET);
            return null;
        }
        else
        {
            double[] newVector = this.vector;
            for (int i = 0; i < this.vectorSize; i++)
                newVector[i] = this.vector[i] + addVector.getVector()[i];
            return new Vector(newVector,this.vectorSize);
        }
    }
    Vector vectorDifference(Vector subtractVector)
    {
        subtractVector = subtractVector.constantMultiplication(-1);
        Vector resultVector;
        resultVector = this.vectorAddition(subtractVector);
        return resultVector;
    }
    Matrix vectorToMatrix()
    {
        double[][] convertMatrix = new double[this.vectorSize][1];
        for (int i = 0; i < this.vectorSize; i++)
            convertMatrix[i][0] = this.vector[i];
        return new Matrix(convertMatrix, this.vectorSize, 1);
    }
}
