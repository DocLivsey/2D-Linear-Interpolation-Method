import java.io.*;
public class Main {
    public static final String RESET = "\u001B[0m";
    public static final String ERROR = "\u001B[31m"; // RED
    public static final String INPUT = "\u001B[32m"; // GREEN
    public static final String COMMENT = "\u001B[33m"; // YELLOW
    public static final String CHOOSE = "\u001B[34m"; // BLUE
    public static final String OUTPUT = "\u001B[35m"; // PURPLE
    public static final String HEADER_OUTPUT = "\u001B[36m"; // CYAN
    public static void main(String[] args) throws FileNotFoundException {
        String pathToVector = "D:/My Files/All Scripts/численные методы/двумерное линейное интерполирование/2D-Linear-Interpolation/java solutions/src/vectorInput.txt";
        String pathToMatrix = "D:/My Files/All Scripts/численные методы/двумерное линейное интерполирование/2D-Linear-Interpolation/java solutions/src/matrixInput.txt";

        Interpolation plot = new Interpolation(pathToVector, pathToVector, pathToMatrix, 1);
        plot.linear2DMethod();

    }
}