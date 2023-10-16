import java.io.*;
import java.util.*;
public class Interpolation {
    protected Vector x;
    protected Vector y;
    protected Matrix values;
    Interpolation(String pathToFileX, String pathToFileY, String pathToFileValues) throws FileNotFoundException {
        this.x = new Vector(pathToFileX);
        this.y = new Vector(pathToFileY);
        this.values = new Matrix(pathToFileValues);
    }
    Interpolation(Vector x, Vector y, Matrix values)
    {
        this.x = x; this.y = y;
        this.values = values;
    }
}
