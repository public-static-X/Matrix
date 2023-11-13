//import java.text.DecimalFormat;
import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

public class Matrix {
    
    //private DecimalFormat df = new DecimalFormat("#.");
    private int nrow;
    private int ncol;
    private double[][] matrix;
    private double[][] Q;

    public Matrix(double[][] matrix) {

        this.nrow = matrix.length;
        this.ncol = matrix[0].length;
        this.matrix = new double[this.nrow][this.ncol];
        
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                this.matrix[i][j] = matrix[i][j];
            }
        }
    }

    public int getNumRows(){
        return this.nrow;
    }

    public int getNumCols(){
        return this.ncol;
    }

    @Override
    public String toString() {
        String printer = "";
        for (int i = 0; i < this.nrow; i++) {
            printer += "|";
            for (int j = 0; j < this.ncol; j++) {
                printer += (this.matrix[i][j] + " ");
            }
        printer += "|\n";
        }
        return printer;
    }


    public Matrix add(Matrix B) {
        double[][] result = new double[this.nrow][this.ncol];

        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                result[i][j] = this.matrix[i][j] + B.matrix[i][j];
            }
        }

        return new Matrix(result);
    }

    public Matrix subtract(Matrix B) {
        double[][] result = new double[this.nrow][this.ncol];

        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                result[i][j] = this.matrix[i][j] - B.matrix[i][j];
            }
        }

        return new Matrix(result);
    }

    public Matrix getRow(int index) {
        double[][] rowMatrix = new double[1][this.ncol];
        for (int i = 0; i < this.ncol; i++) {
            rowMatrix[0][i] = this.matrix[index][i];
        }
        return new Matrix(rowMatrix);
    }

    public Matrix getColumn(int index) {
        double[][] colMatrix = new double[this.nrow][1];
        for (int i = 0; i < this.nrow; i++) {
            colMatrix[i][0] = this.matrix[i][index];
        }
        return new Matrix(colMatrix);
    }

    public Matrix transpose() {
        double[][] temp = new double[this.ncol][this.nrow];

        for (int i = 0; i < this.ncol; i++) {
            for (int j = 0; j < this.nrow; j++) {
                temp[i][j] = this.matrix[j][i];
            }
        }

        return new Matrix(temp);
    }
        
    public static Matrix Identity(int N) {

        double[][] ident = new double[N][N];

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    ident[i][j] = 1;
                } else {
                    ident[i][j] = 0;
                }
            }
        }
        return new Matrix(ident);
    }

    public double getElement(int row, int col) {
        return this.matrix[row][col];
    }

    private double vectorMult(Matrix B) {
        double result = 0;
        for (int i = 0; i < this.ncol; i++) {
            result += this.getElement(0, i) * B.getElement(i, 0);
        }
        return result;
    }

    public double dotProduct(Matrix vector) { // 2 vectors of same array size multiplied.
        double result = 0;

        if(this.isColumnVector() && vector.isColumnVector()) {
            result = this.transpose().vectorMult(vector);
        } else if(this.isRowVector() && vector.isRowVector()) {
            result = this.vectorMult(vector.transpose());
        } else if(this.isRowVector() && vector.isColumnVector()) {
            result = this.vectorMult(vector);
        } else if(this.isColumnVector() && vector.isRowVector()) {
            result = this.transpose().vectorMult(vector.transpose());
        } else {
            throw new IllegalArgumentException("Vectors not conformal for multiplication");
        }

    return result;
    }

    public Matrix multiplyMatrix(Matrix MatrixTwo) {

        double[][] result = new double[this.nrow][MatrixTwo.ncol];

        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < MatrixTwo.ncol; j++) {
                result[i][j] = this.getRow(i).dotProduct(MatrixTwo.getColumn(j));
            }
        }
        return new Matrix(result);
    }

    public Matrix scalarMultiply(double scalar) {
        double[][] result = new double[this.nrow][this.ncol];
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                result[i][j] = this.matrix[i][j] * scalar;
            }
        }
        return new Matrix(result);
    }

    private boolean isRowVector() { // vector needs to be like [1,2,1] for example. must be a row vector
        return this.nrow == 1;
    }

    private boolean isColumnVector() { // vector needs to be like [1,2,1].transpose() or in column format
        return this.ncol == 1;
    }

    

    public double norm() { // maginitude of vector x = ||x||
        double result = 0;
        if(this.isRowVector() || this.isColumnVector()){
        result = this.dotProduct(this);
        } else{
            throw new IllegalArgumentException("Not a vector");
        }
        return Math.sqrt(result);
    }


    public Matrix normalizeVector() { // Mathematically saying x / ||x|| which gives a magnitude of 1
        //double[][] vector = new double[this.nrow][this.ncol];
        if(this.isRowVector() || this.isColumnVector()) {
            return this.scalarMultiply(1 / this.norm()) ;
        } else {
            throw new IllegalArgumentException("Not a vector");
        }
    }

    public static Matrix projection(Matrix u, Matrix a) { // proj_u_(a) = <u,a>/<u,u> .u = scalar multiple of u
        Matrix proj =u.scalarMultiply( (u.dotProduct(a) / u.dotProduct(u)) );
        return proj;
    }

    public Matrix columnBind(Matrix C) {
        
        if(this.nrow != C.nrow) throw new IllegalArgumentException("Incompatible number of rows");

        int numCol = this.ncol + C.ncol;
        double[][] newMat = new double[this.nrow][numCol];
        
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                newMat[i][j] = this.matrix[i][j];
            }
        }

        for (int i = 0; i < this.nrow; i++) {
            for (int j = this.ncol; j < numCol; j++) {
                newMat[i][j] = C.matrix[i][j - this.ncol];
            }
        }

        return new Matrix (newMat);
    }

    public Matrix rowBind(Matrix C) {
        
        if(this.ncol != C.ncol) throw new IllegalArgumentException("Incompatible number of columns");

        int numRow = this.nrow + C.nrow;
        double[][] newMat = new double[numRow][this.ncol];
        
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                newMat[i][j] = this.matrix[i][j];
            }
        }

        for (int i = this.nrow; i < numRow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                newMat[i][j] = C.matrix[i - this.nrow][j];
            }
        }

        return new Matrix (newMat);
    }
    
    public void QRdecomposition() {

        if(this.ncol != this.nrow) throw new IllegalArgumentException("Currently only works for square matrices");

        Matrix U = this.getColumn(0).normalizeVector();
        Matrix column;
        
        for (int i = 1; i < this.ncol; i++) {

            column = this.getColumn(i);
            
            for (int j = 0; j < U.getNumCols(); j++) {
                column = column.subtract( Matrix.projection(U.getColumn(j), column) );
            }
            U = U.columnBind(column.normalizeVector());
        }
        this.Q = U.matrix;
    }

    public Matrix get_Q_matrix() {
        return new Matrix(this.Q);
    }

    public Matrix get_R_matrix() {
        return this.get_Q_matrix().transpose().multiplyMatrix(this);
    }

    public Matrix Round() {
        double[][] roundedMatrix = new double[this.nrow][this.ncol];
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
               roundedMatrix[i][j] = Math.round(100 * this.matrix[i][j] ) / 100;
            }
        }
        return new Matrix(roundedMatrix);
    }

    private Matrix back_substitution() {
        double[][] X = new double[this.nrow][this.ncol];
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.ncol; j++) {
                X[i][j] = 0;
            }
        }

        for (int column = 0; column < this.nrow; column++) {
            double x;
            int counter;

            Matrix I = Matrix.Identity(this.nrow).getColumn(column);
            for (int j = this.nrow - 1; j >= 0; j--) {
                x = I.getElement(j, 0);
                counter = this.nrow - 1;
                while(counter > 0) {
                    x = x - this.matrix[j][counter] * X[counter][column];
                    counter--;
                }
                X[j][column] = x / this.matrix[j][j];
            }
        }
            return new Matrix(X);
    }

    public boolean equals(Matrix B) {

        if(this.nrow != B.nrow | this.ncol != B.ncol) return false;
        for (int i = 0; i < this.nrow; i++) {
            for (int j = 0; j < this.nrow; j++) {
                if(this.matrix[i][j] != B.matrix[i][j]) return false;
            }
        }
        return true;
    }

    public Matrix inverse() {
        this.QRdecomposition();
        Matrix INV =  this.get_R_matrix().back_substitution().multiplyMatrix(get_Q_matrix().transpose());

        if(this.multiplyMatrix(INV).Round().equals(Matrix.Identity(this.nrow))) return INV;
        else throw new IllegalArgumentException("Matrix is singular");
    }

    public static Matrix read_txt(int numRow, int numCol, String filePath) {


        ArrayList<String> arr = new ArrayList<>();
        ArrayList<String> lines = new ArrayList<>();

        //int numRow = 150;
        //int numCol = 6;
        double[][] mat = new double[numRow][numCol];

        File file = new File(filePath);

        try {
            
            int index = 0; // this is row index
            Scanner scan = new Scanner(file);
            
            while(scan.hasNextLine()) {
                arr.add(scan.nextLine());
                Scanner line = new Scanner(arr.get(index));
                while(line.hasNext()) {
                    lines.add( line.next()) ;
                }
                
                for (int i = 0; i < lines.size(); i++) {
                    mat[index][i] =  (Double.parseDouble(lines.get(i)) );
                }
                
                lines.clear();
                index++;
            }
            scan.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    
        Matrix Data = new Matrix(mat);
        return (Data);

    }
    
}
