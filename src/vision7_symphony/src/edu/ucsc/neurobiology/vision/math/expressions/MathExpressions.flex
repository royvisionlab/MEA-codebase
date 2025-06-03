// JFlex specification for a simple expression scanner.
// @author Dumitru Petrusca, University of California, Santa Cruz.

package edu.ucsc.neurobiology.vision.math.expressions;

import java_cup.runtime.*;

%%

%class ExpressionScanner
%unicode
%line
%column
%cup

%{
  StringBuffer string = new StringBuffer();

  public ExpressionScanner(String s) {
    this(new java.io.StringReader(s));
  }

  private Symbol symbol(int type) {
    return new Symbol(type, yyline + 1, yycolumn+1);
  }

  private Symbol symbol(int type, Object value) {
    return new Symbol(type, yyline+1, yycolumn+1, value);
  }
%}


/* main character classes */
LineTerminator = \r|\n|\r\n
WhiteSpace = {LineTerminator} | [ \t\f]

/* identifiers */
Identifier = [:jletter:][:jletterdigit:]*

/* string literals */
StringCharacter = [^\r\n\"\\]

/* floating point literals */
DoubleLiteral = ({FLit1}|{FLit2}|{FLit3}) {Exponent}?

FLit1    = [0-9]+ \. [0-9]*
FLit2    = \. [0-9]+
FLit3    = [0-9]+
Exponent = [eE] [+-]? [0-9]+


%state STRING
%%


<YYINITIAL> {
  /* reserved words */
  "true"                         { return symbol(ExpressionParserSymbols.TRUE); }
  "false"                        { return symbol(ExpressionParserSymbols.FALSE); }
  "pi"                           { return symbol(ExpressionParserSymbols.PI); }
  "e"                            { return symbol(ExpressionParserSymbols.E); }

  /* function names */
  "sin"                          { return symbol(ExpressionParserSymbols.SIN); }
  "cos"                          { return symbol(ExpressionParserSymbols.COS); }
  "exp"                          { return symbol(ExpressionParserSymbols.EXP); }
  "ln"                           { return symbol(ExpressionParserSymbols.LN); }
  "abs"                          { return symbol(ExpressionParserSymbols.ABS); }
  "log10"                        { return symbol(ExpressionParserSymbols.LOG10); }

  "mean"                         { return symbol(ExpressionParserSymbols.MEAN); }
  "sum"                          { return symbol(ExpressionParserSymbols.SUM); }
  "radNorm"                      { return symbol(ExpressionParserSymbols.RAD_NORM); }
  "min"                          { return symbol(ExpressionParserSymbols.MIN); }
  "max"                          { return symbol(ExpressionParserSymbols.MAX); }
  "extreme"                      { return symbol(ExpressionParserSymbols.EXTREME); }
  "norm"             	     	 { return symbol(ExpressionParserSymbols.NORM); }
  "normrms"            	     	 { return symbol(ExpressionParserSymbols.NORMRMS); }
  "fft"                          { return symbol(ExpressionParserSymbols.FFT); }
  "d"                            { return symbol(ExpressionParserSymbols.D); }

  /* separators */
  "("                            { return symbol(ExpressionParserSymbols.LPAREN); }
  ")"                            { return symbol(ExpressionParserSymbols.RPAREN); }
  "["                            { return symbol(ExpressionParserSymbols.LBRACK); }
  "]"                            { return symbol(ExpressionParserSymbols.RBRACK); }
  ","                            { return symbol(ExpressionParserSymbols.COMMA); }

  /* operators */
  "+"                            { return symbol(ExpressionParserSymbols.ADD); }
  "-"                            { return symbol(ExpressionParserSymbols.SUB); }
  "*"                            { return symbol(ExpressionParserSymbols.MUL); }
  "/"                            { return symbol(ExpressionParserSymbols.DIV); }
  "^"                            { return symbol(ExpressionParserSymbols.POW); }
  "%"                            { return symbol(ExpressionParserSymbols.MOD); }
  "#"                            { return symbol(ExpressionParserSymbols.JOIN); }

  "<"                            { return symbol(ExpressionParserSymbols.LESS); }
  "<="                           { return symbol(ExpressionParserSymbols.LESS_EQUAL); }
  ">"                            { return symbol(ExpressionParserSymbols.GREATER); }
  ">="                           { return symbol(ExpressionParserSymbols.GREATER_EQUAL); }
  "!"                            { return symbol(ExpressionParserSymbols.NOT); }
  "&"                            { return symbol(ExpressionParserSymbols.AND); }
  "|"                            { return symbol(ExpressionParserSymbols.OR); }
  "=="                           { return symbol(ExpressionParserSymbols.EQUAL); }
  "!="                           { return symbol(ExpressionParserSymbols.NOT_EQUAL); }

  /* numeric literals */
  {DoubleLiteral}                { return symbol(ExpressionParserSymbols.NUMBER, new Double(yytext())); }

  /* identifiers */
  {Identifier}                   { return symbol(ExpressionParserSymbols.IDENTIFIER, yytext()); }

  /* string */
  \"                             { yybegin(STRING); string.setLength(0); }

  /* whitespace */
  {WhiteSpace}                   { /* ignore */ }
}


<STRING> {
  \"                             { yybegin(YYINITIAL); return symbol(ExpressionParserSymbols.STRING, string.toString()); }
  {StringCharacter}+             { string.append( yytext() ); }

  /* escape sequences */
  "\\b"                          { string.append( '\b' ); }
  "\\t"                          { string.append( '\t' ); }
  "\\n"                          { string.append( '\n' ); }
  "\\f"                          { string.append( '\f' ); }
  "\\r"                          { string.append( '\r' ); }
  "\\\""                         { string.append( '\"' ); }
  "\\'"                          { string.append( '\'' ); }
  "\\\\"                         { string.append( '\\' ); }

  /* error cases */
  \\.                            { throw new RuntimeException("Illegal escape sequence \""+yytext()+"\""); }
  {LineTerminator}               { throw new RuntimeException("Unterminated string at end of line"); }
}


<<EOF>>                          { return symbol(ExpressionParserSymbols.EOF); }
