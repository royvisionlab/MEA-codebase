package edu.ucsc.neurobiology.vision.test;


public class Vector2D {
    public static final int RADIANS = 0;
    public static final int DEGREES = 1;
    private double x, y;


    public Vector2D(double x, double y) {
        this.x = x;
        this.y = y;
    }


    public Vector2D(double phi, double rho, int phiUnits) {
        switch (phiUnits) {
            case RADIANS:
                this.x = rho * Math.cos(phi);
                this.y = rho * Math.sin(phi);
                break;

            case DEGREES:
                this.x = rho * Math.cos(phi * Math.PI / 180.0);
                this.y = rho * Math.sin(phi * Math.PI / 180.0);
                break;

            default:
                throw new IllegalArgumentException("Unknown phi units.");
        }
    }


    public Vector2D add(Vector2D vector) {
        this.x += vector.x;
        this.y += vector.y;

        return this;
    }


    public Vector2D scale(double factor) {
        this.x *= factor;
        this.y *= factor;

        return this;
    }


    public double project(Vector2D vector) {
        return dotProduct(vector) / (Math.sqrt(x * x + y * y));
    }


    public double dotProduct(Vector2D vector) {
        return x * vector.x + y * vector.y;
    }


    public double modulusSquare() {
        return x * x + y * y;
    }


    public double modulus() {
        return Math.sqrt(x * x + y * y);
    }


    public double getPhiRadians() {
        return Math.atan(y / x);
    }


    public double getX() {
        return x;
    }


    public double getY() {
        return y;
    }


    public double getPhiDegrees() {
//        return Math.atan(y/x)*180/Math.PI;
        double phi = Math.atan2(y, x) * 180 / Math.PI;
        while (phi < 0) {
            phi += 360;
        }
        return phi;
    }
}
