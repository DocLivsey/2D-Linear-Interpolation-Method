import java.io.*;
import java.util.*;

public class Interpolation {
    protected ArrayList<Vector> arguments;
    protected Matrix values;
    protected int compactCoefficient;
    protected final double epsilon = 1E-10;
    Interpolation(String pathToFileX, String pathToFileY, String pathToFileValues, int compactCoefficient) throws FileNotFoundException {
        this.arguments = new ArrayList<>();
        this.arguments.add(new Vector(pathToFileX));
        this.arguments.add(new Vector(pathToFileY));
        this.values = new Matrix(pathToFileValues);
        this.compactCoefficient = compactCoefficient;
    }
    Interpolation(Vector x, Vector y, Matrix values, int compactCoefficient)
    {
        this.arguments = new ArrayList<>();
        this.arguments.add(x);
        this.arguments.add(y);
        this.values = values;
        this.compactCoefficient = compactCoefficient;
    }
    ArrayList<Vector> getArguments()
    { return this.arguments; }
    Vector getX()
    { return this.arguments.get(0); }
    Vector getY()
    { return this.arguments.get(1); }
    Matrix getValues()
    { return this.values; }
    double pointStepArithmetic(double step, Vector point)
    { return point.getVector()[0] + step; }
    double methodPointStepArithmetic(Vector point, int expandedSize)
    { return (point.getItem(point.getVectorSize() - 1) - point.getItem(0)) / (expandedSize - 1); }
    void vectorExpansion(Vector points, double newLength, double step)
    {
        Vector partOfPoints;
        for (int i = 0; i < newLength - 1; i++)
        {
            partOfPoints = points.partOfVector(i, i + 1);
            double delta = this.pointStepArithmetic(step, partOfPoints);
            if (!points.isInVector(delta))
                points.addItem(delta); points.sort();
        }
    }
    Matrix matrixExpansion(Matrix valuesMatrix, int coefficient)
    {
        int newRowsCount = valuesMatrix.getRowsCount();
        int newColumnsCount = valuesMatrix.getColumnsCount();
        for (int i = 0; i < coefficient; i++)
        {
            newRowsCount = newRowsCount * 2 - 1;
            newColumnsCount = newColumnsCount * 2 - 1;
        }
        Matrix newValuesMatrix = new Matrix(newRowsCount, newColumnsCount);
        for (int i = 0, placeI = 0; i < valuesMatrix.getRowsCount(); i++, placeI += newRowsCount / (valuesMatrix.getRowsCount() - 1))
            for (int j = 0, placeJ = 0; j < valuesMatrix.getColumnsCount(); j++, placeJ += newColumnsCount / (valuesMatrix.getColumnsCount() - 1))
                newValuesMatrix.setItem(placeI, placeJ, valuesMatrix.getItem(i, j));
        return newValuesMatrix;
    }
    void setLeftUpTriangle(Vector leftUpTriangle, Vector pointsX, Vector pointsY, Matrix partOfValues, Vector xToPlace, Vector yToPlace)
    {
        for (int i = 0, placeInd = 0; i < partOfValues.getRowsCount(); i++)
        {
            for (int j = 0; j < partOfValues.getColumnsCount(); j++)
                if (i != partOfValues.getRowsCount() - 1 || j != partOfValues.getColumnsCount() - 1)
                {
                    leftUpTriangle.setItem(placeInd, partOfValues.getItem(i, j));
                    xToPlace.setItem(placeInd, pointsX.getItem(j));
                    yToPlace.setItem(placeInd, pointsY.getItem(i));
                    placeInd++;
                }
        }
    }
    void setRightDownTriangle(Vector rightDownTriangle, Vector pointsX, Vector pointsY, Matrix partOfValues, Vector xToPlace, Vector yToPlace)
    {
        for (int i = partOfValues.getRowsCount() - 1, placeInd = 0; i >= 0; i--)
        {
            for (int j = partOfValues.getColumnsCount() - 1; j >= 0; j--)
                if (i != 0 || j != 0)
                {
                    rightDownTriangle.setItem(placeInd, partOfValues.getItem(i, j));
                    xToPlace.setItem(placeInd, pointsX.getItem(j));
                    yToPlace.setItem(placeInd, pointsY.getItem(i));
                    placeInd++;
                }
        }
    }
    double planeCoefficient(Vector x, Vector y, Vector z, String coeffName)
    {
        Matrix coeffMatrix = new Matrix(3, 3);
        switch (coeffName)
        {
            case "A":
                coeffMatrix.setColumn(new double[]{1, 1, 1}, 0);
                coeffMatrix.setColumn(y.getVector(), 1);
                coeffMatrix.setColumn(z.getVector(), 2);
                break;
            case "B":
                coeffMatrix.setColumn(x.getVector(), 0);
                coeffMatrix.setColumn(new double[]{1, 1, 1}, 1);
                coeffMatrix.setColumn(z.getVector(), 2);
                break;
            case "C":
                coeffMatrix.setColumn(x.getVector(), 0);
                coeffMatrix.setColumn(y.getVector(), 1);
                coeffMatrix.setColumn(new double[]{1, 1, 1}, 2);
                break;
            case "D":
                coeffMatrix.setColumn(x.getVector(), 0);
                coeffMatrix.setColumn(y.getVector(), 1);
                coeffMatrix.setColumn(z.getVector(), 2);
                break;
        }
        if (Math.abs(coeffMatrix.matrix3By3Determinant()) < this.epsilon)
        { return 0; }
        else
        {
            if (coeffName.equals("D"))
                return -coeffMatrix.matrix3By3Determinant();
            return coeffMatrix.matrix3By3Determinant();
        }
    }
    void valueArithmetic(Vector targetX, double currentY, int rowPlaceInd, int colPlaceInd, Matrix changeValue, double detB, double detA, double detD, double detC, int placeI, int placeJ, int rowsSectionLength, int colsSectionLength)
    {
        double currentX;
        double value;
        currentX = targetX.getItem(placeJ);
        value = - ((detA * currentX + detB * currentY + detD) / detC);
        if (Double.isNaN(changeValue.getItem(placeI, placeJ)))
        {
            //System.out.println("i= " + placeI + " j= " + placeJ + " borders: " + rowsSectionLength + " " + colsSectionLength);
            changeValue.setItem(placeI, placeJ, value);
        }
    }
    void triangleArithmetic(Vector originX, Vector originY, Matrix partOfValues, Vector targetX, Vector targetY, Matrix changeValue, int rowPlaceInd, int colPlaceInd, int rowsSectionLength, int colsSectionLength)
    {
        Vector xToPlace = new Vector(3);
        Vector yToPlace = new Vector(3);
        Vector leftUpTriangle = new Vector(3);
        this.setLeftUpTriangle(leftUpTriangle, originX, originY, partOfValues, xToPlace, yToPlace);

        double detA = this.planeCoefficient(xToPlace, yToPlace, leftUpTriangle, "A");
        double detB = this.planeCoefficient(xToPlace, yToPlace, leftUpTriangle, "B");
        double detC = this.planeCoefficient(xToPlace, yToPlace, leftUpTriangle, "C");
        double detD = this.planeCoefficient(xToPlace, yToPlace, leftUpTriangle, "D");

        //System.out.println("nodeCoordinates: " + rowPlaceInd + " " + colPlaceInd);

        double currentY;
        for (int placeI = rowPlaceInd; placeI < rowsSectionLength + rowPlaceInd; placeI++)
        {
            currentY = targetY.getItem(placeI);
            for (int placeJ = colPlaceInd; placeJ < colsSectionLength + colPlaceInd - placeI + rowPlaceInd; placeJ++)
            { this.valueArithmetic(targetX, currentY, rowPlaceInd, colPlaceInd, changeValue, detB, detA, detD, detC, placeI, placeJ, rowsSectionLength, colsSectionLength); }

        }

        Vector rightDownTriangle = new Vector(3);
        this.setRightDownTriangle(rightDownTriangle, originX, originY, partOfValues, xToPlace, yToPlace);

        detA = this.planeCoefficient(xToPlace, yToPlace, rightDownTriangle, "A");
        detB = this.planeCoefficient(xToPlace, yToPlace, rightDownTriangle, "B");
        detC = this.planeCoefficient(xToPlace, yToPlace, rightDownTriangle, "C");
        detD = this.planeCoefficient(xToPlace, yToPlace, rightDownTriangle, "D");

        for (int placeI = rowPlaceInd; placeI < rowsSectionLength + rowPlaceInd + 1; placeI++)
        {
            currentY = targetY.getItem(placeI);
            for (int placeJ = colsSectionLength + colPlaceInd; placeJ >= colsSectionLength + colPlaceInd - placeI + rowPlaceInd; placeJ--)
            { this.valueArithmetic(targetX, currentY, rowPlaceInd, colPlaceInd, changeValue, detB, detA, detD, detC, placeI, placeJ, rowsSectionLength, colsSectionLength); }
        }
    }
    void linear2DMethod() throws IOException
    {
        // копируем векторы Х и У, и матрицу значений
        Vector copyX = this.getX().cloneVector();
        Vector copyY = this.getY().cloneVector();
        Matrix copyValues = this.matrixExpansion(this.values, this.compactCoefficient);

        // разделяем вектора значений на части по два на каждом шаге записывая в список
        Vector xPart;
        Vector yPart;

        // шаг интерполирования с учетом коэффициента уплотнения
        double stepY = this.methodPointStepArithmetic(copyY, copyValues.getRowsCount());
        double stepX = this.methodPointStepArithmetic(copyX, copyValues.getColumnsCount());
        this.vectorExpansion(copyY, copyValues.getRowsCount(), stepY);
        this.vectorExpansion(copyX, copyValues.getColumnsCount(), stepX);

        copyX.printVector(); copyX.printVector();
        int rowsSectionLength = copyValues.getRowsCount() / (this.values.getRowsCount() - 1);
        int colsSectionLength = copyValues.getColumnsCount() / (this.values.getColumnsCount() - 1);
        for (int rowInd = 0, newRowInd = 0; rowInd < this.values.getRowsCount() - 1; rowInd++, newRowInd += rowsSectionLength)
        {
            yPart = this.getY().partOfVector(rowInd, rowInd + 1);
            for (int colInd = 0, newColInd = 0; colInd < this.values.getColumnsCount() - 1; colInd++, newColInd += colsSectionLength)
            {
                xPart = this.getX().partOfVector(colInd, colInd + 1);
                Matrix valuesPart = this.values.partOfMatrix(colInd, colInd + 1, rowInd, rowInd + 1);
                this.triangleArithmetic(xPart, yPart, valuesPart, copyX, copyY, copyValues, newRowInd, newColInd, rowsSectionLength, colsSectionLength);
            }
        }
        this.arguments.set(0, copyX); this.arguments.set(1, copyY); this.values = copyValues.cloneMatrix();
        this.values.writeInFile("src/matrixOutput.txt");
        this.arguments.get(0).writeInFile("src/vectorXOutput.txt");
        this.arguments.get(1).writeFormattedInFile("src/vectorYOutput.txt");
        this.values.printFormattedMatrix();
    }

    void advancedLinear2DMethod() throws IOException
    {
        Vector copyX = this.getX().cloneVector();
        Vector copyY = this.getY().cloneVector();
        Matrix copyValues = this.matrixExpansion(this.values, 1);

        for (int it = 0; it < this.compactCoefficient; it++)
        {
            Vector xPart;
            Vector yPart;
            double dX = this.methodPointStepArithmetic(copyX, copyValues.getColumnsCount());
            double dY = this.methodPointStepArithmetic(copyY, copyValues.getRowsCount());
            this.vectorExpansion(copyX, copyValues.getColumnsCount(), dX);
            this.vectorExpansion(copyY, copyValues.getRowsCount(), dY);


            this.values = copyValues.cloneMatrix();
            this.arguments.set(0, copyX);
            this.arguments.set(1, copyY);

            copyValues = this.matrixExpansion(copyValues, 1);
        }
        this.values.printMatrix();
    }
}
