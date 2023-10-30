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
    double methodStepArithmetic(Vector point, int expandedSize)
    {
        int coefficient = expandedSize / 2;
        return (point.getItem(1) - point.getItem(0)) / coefficient;
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
        if (Math.abs(coeffMatrix.gaussianDeterminant()) < this.epsilon)
        { return 0; }
        else
        {
            if (coeffName.equals("D"))
                return -coeffMatrix.gaussianDeterminant();
            return coeffMatrix.gaussianDeterminant();
        }
    }
    void triangleValueArithmetic(Vector originX, Vector originY, Matrix partOfValues, Vector targetX, Vector targetY, Matrix changeValue, int rowPlaceInd, int colPlaceInd)
    {
        Vector xToPlace = new Vector(3);
        Vector yToPlace = new Vector(3);
        Vector leftUpTriangle = new Vector(3);
        this.setLeftUpTriangle(leftUpTriangle, originX, originY, partOfValues, xToPlace, yToPlace);
        leftUpTriangle.printFormattedVector(); xToPlace.printFormattedVector(); yToPlace.printFormattedVector();

        double detA = this.planeCoefficient(xToPlace, yToPlace, leftUpTriangle, "A");
        double detB = this.planeCoefficient(xToPlace, yToPlace, leftUpTriangle, "B");
        double detC = this.planeCoefficient(xToPlace, yToPlace, leftUpTriangle, "C");
        double detD = this.planeCoefficient(xToPlace, yToPlace, leftUpTriangle, "D");

        System.out.println(rowPlaceInd + " " + colPlaceInd);
        double currentY, currentX, value;
        for (int placeI = rowPlaceInd; placeI < changeValue.getRowsCount()/2 + rowPlaceInd; placeI++)
        {
            currentY = targetY.getItem(placeI);
            for (int placeJ = colPlaceInd; placeJ < changeValue.getColumnsCount()/2 + colPlaceInd - placeI + rowPlaceInd + 1; placeJ++)
            {
                currentX = targetX.getItem(placeJ);
                value = - ((detA * currentX + detB * currentY + detD) / detC);
                if (changeValue.getItem(placeI, placeJ) == 0)
                    changeValue.setItem(placeI, placeJ, value);
            }
        }
        changeValue.printMatrix();

        Vector rightDownTriangle = new Vector(3);
        this.setRightDownTriangle(rightDownTriangle, originX, originY, partOfValues, xToPlace, yToPlace);
        rightDownTriangle.printFormattedVector(); xToPlace.printFormattedVector(); yToPlace.printFormattedVector();

        detA = this.planeCoefficient(xToPlace, yToPlace, rightDownTriangle, "A");
        detB = this.planeCoefficient(xToPlace, yToPlace, rightDownTriangle, "B");
        detC = this.planeCoefficient(xToPlace, yToPlace, rightDownTriangle, "C");
        detD = this.planeCoefficient(xToPlace, yToPlace, rightDownTriangle, "D");

        Matrix changeMatrixPart = changeValue.partOfMatrix(colPlaceInd, colPlaceInd + changeValue.getColumnsCount()/2, rowPlaceInd, rowPlaceInd + changeValue.getRowsCount()/2);
        changeMatrixPart.printMatrix();

        for (int placeI = rowPlaceInd + 1; placeI < changeValue.getRowsCount()/2 + rowPlaceInd + 1; placeI++)
        {
            currentY = targetY.getItem(placeI);
            for (int placeJ = changeValue.getColumnsCount() / 2 + colPlaceInd; placeJ >= changeValue.getColumnsCount() / 2 + colPlaceInd - placeI + rowPlaceInd; placeJ--)
            {
                currentX = targetX.getItem(placeJ);
                value = - ((detA * currentX + detB * currentY + detD) / detC);
                if (changeValue.getItem(placeI, placeJ) == 0)
                    changeValue.setItem(placeI, placeJ, value);
            }
        }
        changeValue.printMatrix();
    }
    void vectorExpansion(Vector points, Vector partOfPoints, double step)
    {
        double delta = this.pointStepArithmetic(step, partOfPoints);
        if (!points.isInVector(delta))
            points.addItem(delta); points.sort();
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
        for (int i = 0, placeI = 0; placeI < newRowsCount; i++, placeI+=(newRowsCount - valuesMatrix.getRowsCount()) / 2 + 1)
            for (int j = 0, placeJ = 0; placeJ < newColumnsCount; j++, placeJ+=(newColumnsCount - valuesMatrix.getColumnsCount()) / 2 + 1)
                newValuesMatrix.setItem(placeI, placeJ, valuesMatrix.getItem(i, j));
        return newValuesMatrix;
    }
    void linear2DMethod()
    {
        // начальные значения границ интерполирования
        double left = this.getX().getItem(0);
        double right = this.getY().getItem(this.getX().getVectorSize() - 1);

        // копируем векторы Х и У, и матрицу значений
        Vector copyX = this.getX().cloneVector();
        Vector copyY = this.getY().cloneVector();
        Matrix copyValues = this.matrixExpansion(this.values, this.compactCoefficient);

        // разделяем вектора значений на части по два на каждом шаге записывая в список
        Vector xPart;
        Vector yPart;

        // шаг интерполирования с учетом коэффициента уплотнения
        double stepY = this.methodStepArithmetic(copyY, copyValues.getRowsCount());
        double stepX = this.methodStepArithmetic(copyX, copyValues.getColumnsCount());
        for (int i = 0; i < copyValues.getRowsCount() - 1; i++)
        {
            yPart = copyY.partOfVector(i, i + 1);
            this.vectorExpansion(copyY, yPart, stepY);
            for (int j = 0; j < copyValues.getColumnsCount() - 1; j++)
            {
                xPart = copyX.partOfVector(j, j + 1);
                this.vectorExpansion(copyX, xPart, stepX);
            }
        }
        for (int rowInd = 0, newRowInd = 0; rowInd < this.values.getRowsCount() - 1; rowInd++, newRowInd+=copyValues.getRowsCount()/2)
        {
            yPart = this.getY().partOfVector(rowInd, rowInd + 1);
            for (int colInd = 0, newColInd = 0; colInd < this.values.getColumnsCount() - 1; colInd++, newColInd+=copyValues.getColumnsCount()/2)
            {
                xPart = this.getX().partOfVector(colInd, colInd + 1);
                Matrix valuesPart = this.values.partOfMatrix(colInd, colInd + 1, rowInd, rowInd + 1);

                this.triangleValueArithmetic(xPart, yPart, valuesPart, copyX, copyY, copyValues, newRowInd, newColInd);
            }
        }
        copyY.printVector(); copyX.printVector();
    }
}
