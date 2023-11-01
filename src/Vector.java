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
    double getItem(int index) { return this.vector[index]; }
    void setItem(int index, double replaceItem)
    { this.vector[index] = replaceItem; }
    void setVector(double[] vector)
    {
        System.out.println(Main.CHOOSE + "Вы уверены, что хотите заменить вектор?" + Main.RESET);
        this.vector = vector;
    }
    public boolean equals(Object obj) {
        return super.equals(obj);
    }
    Vector cloneVector()
    { return new Vector(this.vector, this.vectorSize); }
    void printVector()
    {
        System.out.print(Main.HEADER_OUTPUT + "\nВектор размерностью " + vectorSize + Main.OUTPUT + ": \n { ");
        for (int i = 0; i < vectorSize; i++)
        {
            System.out.print(this.getItem(i) + "; ");
        }
        System.out.println("}" + Main.RESET);
    }
    void printFormattedVector()
    {
        System.out.print(Main.HEADER_OUTPUT + "\nВектор размерностью " + vectorSize + Main.OUTPUT + ": \n { ");
        for (int i = 0; i < vectorSize; i++)
        {
            DecimalFormat shortOut = new DecimalFormat("#.##");
            String result = shortOut.format(this.getItem(i));
            System.out.print(result + "; ");
        }
        System.out.println("}" + Main.RESET);
    }
    void writeInFile(String pathToFile) throws IOException {
        FileWriter fileWriter = new FileWriter(pathToFile);
        for (int i = 0; i < this.vectorSize; i++)
            fileWriter.write(this.getItem(i) + " ");
        fileWriter.close();
    }
    void writeFormattedInFile(String pathToFile) throws IOException {
        FileWriter fileWriter = new FileWriter(pathToFile);
        for (int i = 0; i < this.vectorSize; i++)
        {
            DecimalFormat shortOut = new DecimalFormat("#.##");
            String result = shortOut.format(this.getItem(i));
            fileWriter.write(result + " ");
        }
        fileWriter.close();
    }
    void addItem(double item)
    {
        this.vectorSize ++;
        double[] newVector = new double[this.vectorSize];
        for (int i = 0; i < this.vectorSize; i++)
        {
            if (i != this.vectorSize - 1)
                newVector[i] = this.getItem(i);
            else
                newVector[i] = item;
        }
        this.vector = newVector;
    }
    void addItemBefore(double item, int index)
    {
        this.vectorSize ++;
        double[] newVector = new double[this.vectorSize];
        newVector[index] = item;
        for (int i = 0; i < this.vectorSize; i++)
        {
            if (i < index)
                newVector[i] = this.getItem(i);
            else if (i > index)
                newVector[i] = this.getItem(i - 1);
        }
        this.vector = newVector;
    }
    void addItemAfter(double item, int index)
    {
        this.vectorSize ++;
        double[] newVector = new double[this.vectorSize];
        newVector[index + 1] = item;
        for (int i = 0; i < this.vectorSize; i++)
        {
            if (i < index + 1)
                newVector[i] = this.getItem(i);
            else if (i > index + 1)
                newVector[i] = this.getItem(i - 1);
        }
        this.vector = newVector;
    }
    Vector partOfVector(int leftBorder, int rightBorder)
    {
        int newVectorSize = rightBorder - leftBorder + 1;
        double[] vectorPart = new double[newVectorSize];
        for (int oldIndex = leftBorder, newIndex = 0; oldIndex < rightBorder + 1; oldIndex++, newIndex++)
            vectorPart[newIndex] = this.getItem(oldIndex);
        return new Vector(vectorPart, newVectorSize);
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
                newVector[i] = this.getItem(i) + addVector.getItem(i);
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
            convertMatrix[i][0] = this.getItem(i);
        return new Matrix(convertMatrix, this.vectorSize, 1);
    }
    boolean isInVector(double item)
    {
        for (double i : this.vector)
            if (i == item) return true;
        return false;
    }
    void sort()
    { Arrays.sort(this.vector); }
}
