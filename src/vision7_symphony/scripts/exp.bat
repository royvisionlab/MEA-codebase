@echo off
cls
cd src\edu\ucsc\neurobiology\vision\math\expressions
java -cp f:\java-library\java_cup.jar java_cup.Main -parser ExpressionParser -symbols ExpressionParserSymbols < MathExpressions.cup
java -Xmx128m -jar f:\java-library\JFlex\lib\JFlex.jar MathExpressions.flex
cd ..\..\..\..\..\..\..