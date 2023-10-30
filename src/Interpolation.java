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
    void triangleValueArithmetic(Vector originX, Vector originY, Matrix partOfValues, Vector targetX, Vector targetY, Matrix changeValue, int rowPlaceInd, int colPlaceInd)
    {
        Vector lefUpTriangle = new Vector(3);
        Vector rightUpTriangle = new Vector(3);
        // верхний левый треугольник
        for (int i = 0, placeInd = 0; i < partOfValues.getRowsCount(); i++)
        {
            for (int j = partOfValues.getColumnsCount() - 1; j >= 0; j--)
                if (i != partOfValues.getRowsCount() - 1 || j != partOfValues.getColumnsCount() - 1)
                {
                    lefUpTriangle.setItem(placeInd, partOfValues.getItem(i, j));
                    placeInd++;
                }
        }
        // нижний правый треугольник
        for (int i = partOfValues.getRowsCount() - 1, placeInd = 0; i >= 0; i--)
        {
            for (int j = 0; j < partOfValues.getColumnsCount(); j++)
                if (i != 0 || j != 0)
                {
                    rightUpTriangle.setItem(placeInd, partOfValues.getItem(i, j));
                    placeInd++;
                }
        }
        lefUpTriangle.printVector();

        Matrix D0 = new Matrix(3, 3);
        D0.setColumn(lefUpTriangle.getVector(), 2);
        for (int i = 0, placeInd = 0; i < originY.getVectorSize(); i++)
        {
            for (int j = originX.getVectorSize() - 1; j >= 0; j--)
                if (i != partOfValues.getRowsCount() - 1 || j != partOfValues.getColumnsCount() - 1)
                {
                    D0.setItem(placeInd, 0, originX.getItem(j));
                    D0.setItem(placeInd, 1, originY.getItem(i));
                    placeInd++;
                }
        }
        double detD0 = D0.gaussianDeterminant();
        Matrix D1 = new Matrix(3, 3);
        D1.setColumn(lefUpTriangle.getVector(), 1);
        for (int i = 0, placeInd = 0; i < originY.getVectorSize(); i++)
        {
            for (int j = originX.getVectorSize() - 1; j >= 0; j--)
                if (i != partOfValues.getRowsCount() - 1 || j != partOfValues.getColumnsCount() - 1)
                {
                    D1.setItem(placeInd, 2, 1);
                    D1.setItem(placeInd, 0, originY.getItem(i));
                    placeInd++;
                }
        }
        double detD1 = D1.gaussianDeterminant();
        Matrix D2 = new Matrix(3, 3);
        D2.setColumn(lefUpTriangle.getVector(), 1);
        for (int i = 0, placeInd = 0; i < originY.getVectorSize(); i++)
        {
            for (int j = originX.getVectorSize() - 1; j >= 0; j--)
                if (i != partOfValues.getRowsCount() - 1 || j != partOfValues.getColumnsCount() - 1)
                {
                    D2.setItem(placeInd, 0, originX.getItem(j));
                    D2.setItem(placeInd, 2, 1);
                    placeInd++;
                }
        }
        double detD2 = D2.gaussianDeterminant();
        Matrix D3 = new Matrix(3, 3);
        for (int i = 0, placeInd = 0; i < originY.getVectorSize(); i++)
        {
            for (int j = originX.getVectorSize() - 1; j >= 0; j--)
                if (i != partOfValues.getRowsCount() - 1 || j != partOfValues.getColumnsCount() - 1)
                {
                    D3.setItem(placeInd, 0, originX.getItem(j));
                    D3.setItem(placeInd, 1, originY.getItem(i));
                    D3.setItem(placeInd, 2, 1);
                    placeInd++;
                }
        }
        double detD3 = D3.gaussianDeterminant();

        System.out.println(rowPlaceInd + " " + colPlaceInd);
        double currentY, currentX, value;
        for (int placeI = rowPlaceInd; placeI < changeValue.getRowsCount()/2 + rowPlaceInd; placeI++)
        {
            currentY = targetY.getItem(placeI);
            for (int placeJ = colPlaceInd; placeJ < changeValue.getColumnsCount()/2 + colPlaceInd - placeI + rowPlaceInd + 1; placeJ++)
            {
                currentX = targetX.getItem(placeJ);
                value = (detD0 - detD1 * currentX - detD2 * currentY) / detD3;
                if (changeValue.getItem(placeI, placeJ) == 0)
                    changeValue.setItem(placeI, placeJ, value);
            }
        }
        changeValue.printMatrix();

        D0.setColumn(lefUpTriangle.getVector(), 2);
        for (int i = originY.getVectorSize() - 1, placeInd = 0; i >= 0; i--)
        {
            for (int j = 0; j < originX.getVectorSize(); j++)
                if (i != 0 || j != 0)
                {
                    D0.setItem(placeInd, 0, originX.getItem(j));
                    D0.setItem(placeInd, 1, originY.getItem(i));
                    placeInd++;
                }
        }
        detD0 = D0.gaussianDeterminant();
        D1.setColumn(lefUpTriangle.getVector(), 1);
        for (int i = originY.getVectorSize() - 1, placeInd = 0; i >= 0; i--)
        {
            for (int j = 0; j < originX.getVectorSize(); j++)
                if (i != 0 || j != 0)
                {
                    D1.setItem(placeInd, 2, 1);
                    D1.setItem(placeInd, 0, originY.getItem(i));
                    placeInd++;
                }
        }
        detD1 = D1.gaussianDeterminant();
        D2.setColumn(lefUpTriangle.getVector(), 1);
        for (int i = originY.getVectorSize() - 1, placeInd = 0; i >= 0; i--)
        {
            for (int j = 0; j < originX.getVectorSize(); j++)
                if (i != 0 || j != 0)
                {
                    D2.setItem(placeInd, 2, 1);
                    D2.setItem(placeInd, 0, originX.getItem(j));
                    placeInd++;
                }
        }
        detD2 = D2.gaussianDeterminant();
        for (int i = originY.getVectorSize() - 1, placeInd = 0; i >= 0; i--)
        {
            for (int j = 0; j < originX.getVectorSize(); j++)
                if (i != 0 || j != 0)
                {
                    D3.setItem(placeInd, 2, 1);
                    D3.setItem(placeInd, 0, originY.getItem(i));
                    D3.setItem(placeInd, 0, originX.getItem(j));
                    placeInd++;
                }
        }
        detD3 = D3.gaussianDeterminant();

        Matrix changeMatrixPart = changeValue.partOfMatrix(colPlaceInd, colPlaceInd + changeValue.getColumnsCount()/2, rowPlaceInd, rowPlaceInd + changeValue.getRowsCount()/2);
        changeMatrixPart.printMatrix();
        for (int i = 1; i < changeMatrixPart.getRowsCount(); i++)
            for (int j = changeMatrixPart.getColumnsCount() - 1; j >= changeMatrixPart.getColumnsCount() - i; j--)
                break;

        for (int placeI = rowPlaceInd + 1; placeI < changeValue.getRowsCount()/2 + rowPlaceInd + 1; placeI++)
        {
            currentY = targetY.getItem(placeI);
            for (int placeJ = changeValue.getColumnsCount() / 2 + colPlaceInd; placeJ >= changeValue.getColumnsCount() / 2 + colPlaceInd - placeI + rowPlaceInd; placeJ--)
            {
                currentX = targetX.getItem(placeJ);
                value = (detD0 - detD1 * currentX - detD2 * currentY) / detD3;
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


    double intermediateValuesArithmetic(double z1, double z2, double pointStep, Vector point)
    {
        double firstMultiplier = (point.getVector()[1] - pointStep) / (point.getVector()[1] - point.getVector()[0]);
        double secondMultiplier = (pointStep - point.getVector()[0]) / (point.getVector()[1] - point.getVector()[0]);

        return firstMultiplier * z1 + secondMultiplier * z2;
    }

    void biLinearMethod()
    {
        // начальные значения шага и границ интерполирования
        double left = this.getX().getVector()[0];
        double right = this.getX().getVector()[this.getX().getVectorSize() - 1];

        for (int it = 0; it < this.compactCoefficient; it++) {
            // копируем векторы Х и У, и матрицу значений
            Vector copyX = this.getX().cloneVector();
            Vector copyY = this.getY().cloneVector();
            Matrix copyValues = this.values.cloneMatrix();

            double stepX = (copyX.getItem(1) - copyX.getItem(0)) / 2;
            double stepY = (copyY.getItem(1) - copyY.getItem(0)) / 2;

            // разделяем вектора значений на части по два на каждом шаге записывая в список
            Vector xToChange;
            Vector yToChange;
            for (int i = 0, placeToAddX = 0; i < this.values.getRowsCount() - 1; i++, placeToAddX += 2) {
                // создаем временный столбец новых промежуточных значений для добавления в исходную матрицу зачений
                double[] tmpCol = new double[this.values.getRowsCount()];
                // новое промежуточное значение Х
                double targetX;
                for (int j = 0; j < this.values.getColumnsCount() - 1; j++) {
                    xToChange = this.getX().partOfVector(j, j + 1);
                    // задаем новые промежуточные значения Х учитывая шаг
                    targetX = this.pointStepArithmetic(stepX, xToChange);
                    // разделяем матрицу значений на части для вычисления промежуточных значений
                    Matrix valuesToChange;
                    if (j == this.values.getColumnsCount() - 2)
                        valuesToChange = this.values.partOfMatrix(i, i + 1, j + 1, j + 1);
                    else
                        valuesToChange = this.values.partOfMatrix(i, i + 1, j, j + 1);

                    // рассчет промежуточных значений
                    double z11, z21, z12, z22, interNodeValueByX1 = 0, interNodeValueByX2 = 0;
                    if (valuesToChange.getRowsCount() == 1) {
                        z11 = valuesToChange.getItem(0, 0);
                        z21 = valuesToChange.getItem(0, 1);
                        interNodeValueByX1 = this.intermediateValuesArithmetic(z11, z21, targetX, xToChange);

                        tmpCol[j + 1] = interNodeValueByX1;
                    } else if (valuesToChange.getRowsCount() == 2) {
                        z11 = valuesToChange.getItem(0, 0);
                        z21 = valuesToChange.getItem(0, 1);
                        interNodeValueByX1 = this.intermediateValuesArithmetic(z11, z21, targetX, xToChange);

                        z12 = valuesToChange.getItem(1, 0);
                        z22 = valuesToChange.getItem(1, 1);
                        interNodeValueByX2 = this.intermediateValuesArithmetic(z12, z22, targetX, xToChange);

                        tmpCol[j] = interNodeValueByX1;
                        tmpCol[j + 1] = interNodeValueByX2;
                    }
                    if (!copyX.isInVector(targetX))
                        copyX.addItem(targetX);
                    copyX.sort();
                }
                copyValues.addColumnAfter(tmpCol, placeToAddX);
            }
            this.arguments.set(0, copyX.cloneVector());
            this.values = copyValues.cloneMatrix();

            // новое промежуточное значение У
            double targetY;
            // создаем временную строку новых промежуточных значений для добавления в исходную матрицу зачений
            double[] tmpRow = new double[this.values.getColumnsCount()];
            for (int i = 0, placeToAddY = 0; i < this.values.getRowsCount() - 1; i++, placeToAddY += 2) {
                yToChange = this.getY().partOfVector(i, i + 1);
                // задаем новые промежуточные значения У учитывая шаг
                targetY = this.pointStepArithmetic(stepY, yToChange);
                Matrix valuesToChange = this.values.partOfMatrix(0, this.values.getColumnsCount() - 1, i, i + 1);
                for (int j = 0; j < valuesToChange.getColumnsCount(); j++) {
                    tmpRow[j] = intermediateValuesArithmetic(valuesToChange.getItem(0, j), valuesToChange.getItem(1, j), targetY, yToChange);
                }
                if (!copyY.isInVector(targetY))
                    copyY.addItemAfter(targetY, placeToAddY);
                copyY.sort();
                copyValues.addRowAfter(tmpRow, placeToAddY);
            }
            copyX.printVector();
            copyY.printVector();
            copyValues.printMatrix();
            this.arguments.set(0, copyX);
            this.arguments.set(1, copyY);
            this.values = copyValues;
        }
    }
}
