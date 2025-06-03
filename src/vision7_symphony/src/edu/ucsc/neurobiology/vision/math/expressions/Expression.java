package edu.ucsc.neurobiology.vision.math.expressions;

import java.util.*;


/**
 * Implements a mathematical expression in the Vision embedded language. The expression
 * is dynamically typed, can have variables and can be executed. For a description of
 * the language look at MathExpressions.flex and MathExpressions.cup.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Expression {
    private String expression;


    public Expression(String expression) {
        this.expression = expression;
    }


    /**
     * Evaluate the expression given the name=value variables. Can be called more that once.
     *
     * @param variables HashMap
     * @return Object
     * @throws Exception
     */
    public Object evaluate(HashMap variables) throws Exception {
        ExpressionScanner scanner = new ExpressionScanner(expression);
        ExpressionParser p = new ExpressionParser(scanner, variables);
        try {
            p.parse();
        } catch (NullPointerException e) {
        }
        return p.result;
    }

}
