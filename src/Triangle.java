import java.text.*;
public class Triangle {
    protected double a, b, c;
    Triangle (double a, double b, double c)
    {
        this.a = a;
        this.b = b;
        this.c = c;
    }
    double getA() { return this.a; }
    double getB() { return this.b; }
    double getC() { return this.c; }
    void setA(double a) { this.a = a; }
    void setB(double a) { this.a = a; }
    void setC(double a) { this.a = a; }
    public boolean equals(Object obj)
    { return super.equals(obj); }
    Triangle cloneTriangle()
    { return new Triangle(this.a, this.b, this.c); }
    public void printTriangle()
    {
        System.out.println(Main.HEADER_OUTPUT + "Треугольник со сторонами:" + Main.RESET);
        System.out.println(Main.OUTPUT + "{ " + a + "; " + b + "; " + c + " }" + Main.RESET);
    }
    public void printFormattedTriangle()
    {
        System.out.println(Main.HEADER_OUTPUT + "Треугольник со сторонами:" + Main.RESET);
        DecimalFormat formattedOut = new DecimalFormat("#.##");
    }
}
