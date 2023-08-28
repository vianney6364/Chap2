#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
DataFrame dyn_dispatch_prod_rect_morta_rcpp(DataFrame ET_T1, DataFrame ET_T2) {
	// ici on définit une fonction qui s'appelle myFunction,
	//  qui renvoie un DataFrame (type d'objet de C++/Rcpp pratique pour les interactions avec R)
	//  c'est pour ça qu'on écrit "DataFrame myFunction(...)"
	//  Les arguments sont ensuite ET_T1 et ET_T2 qui sont deux DataFrame
	
	NumericVector T1_inf = ET_T1["TLinf"];
	NumericVector T2_inf = ET_T2["TLinf"];
	NumericVector T1_sup = ET_T1["TLsup"];
	NumericVector T2_sup = ET_T2["TLsup"];
	NumericVector T1_prod = ET_T1["production"];
	// NumericVector c'est un type de C++/Rcpp qui correspond à un vector de nombres réels
	// pour un vecteur d'entier il faudrait écrire IntegerVector T1_x
	// Ici on crée un pointeur (pointer en anglais) vers la colonne "x" du DataFrame "ET_T1", 
	// pour plus de praticité dans la manipulation mais on pourrait directement manipuler le DataFrame
	// Attention c'est un pointeur donc si tu modifie T1_x, ET_T1["x"] sera modifié aussi.
	// L'avantage c'est que ça ne coûte presque rien en terme de vitesse de calcul est de mémoire.
	// 
	// Si on voulait copier la colonne "x" dans un nouvel objet 
	// que l'on peut manipuler indépendemment il faudrait écrire:
	// NumericVector T1_x = clone(ET_T1["x"]);
	// ce qui correspond au comportement de R du type "my_x <- df$x"
	// 
	// À Noter : à la fin de chaque instruction C++ il faut mettre un ";"
	// 
	
	// On créé des raccourcis (pointeurs) pour chacune des autres colonnes
	
	NumericVector T2_prod(ET_T2.nrows()); 
	T2_prod[0]=R_NaN;
	int j=1;
	// Ci-dessus on créé un vecteur du nom de my_prod et de longueur ET_T2.nrows()
	// pour obtenir le nombre de ligne de ET_T2 on écrit ET_T2.nrows() plutôt que nrows(ET_T2)
	// on pourrait aussi stocker directement le nombre de ligne de ET_T2 dans une variable:
	// double matrix_length = ET_T2.nrows();
	// pour éviter de répéter cette opération
	// "double" est un type classique du C++ qui correspond à un réel.
	// Rcout << "je suis pas encore dans la boucle";
	
	for(int i = 1; i < T2_sup.size(); ++i) {
		// syntaxe de la boucle 'for' avec trois éléments
		// int i = 0; définit l'entier sur lequel on boucle
		// i <T2_x.size(); définit la condition de continuité de la boucle. 
		// Ici on continue tant que i est plus petit que la taille du vecteur T2_x
		// ++i définit ce qu'il se passe à chaque itération de la boucle.
		// ici on incrémente juste i de 1
		// 
		// au passage ++i ou i++ correspond à:
		// i = i+1 
		// mais en raccourci
// Rcout << "je suis à tl=" << i << "\n";
		
		// ***** NOTE IMPORTANTE *****
		// les indices commence à 0 en C++ et à 1 en R
		// donc en R x[1] c'est équivalent en C++ à x[0]
		while (T2_sup[i] >= T1_sup[j]) {
		  T2_prod[i]=T2_prod[i]+(T1_sup[j]-max(NumericVector::create(T1_inf[j], T2_inf[i])))/
		    (T1_sup[j]-T1_inf[j])*T1_prod[j];
		  // Rcout << "The value of b : " << T2_prod[i] << "\n";
		  
		  j++;
		  if(j>(T1_inf.size()-1)){
		    break;
		  }
		}
		if (j<=(T1_inf.size()-1)){
		  T2_prod[i]=T2_prod[i]+(T2_sup[i]-max(NumericVector::create(T1_inf[j], T2_inf[i])))/
		    (T1_sup[j]-T1_inf[j])*T1_prod[j];
		}
	}
	// 
	// 	my_prod[i] += T1_x[i]+T1_y[i];
	// 	// my_prod += 1; est équivalent à :
	// 	// my_prod = my_prod +1;
	// 	// c'est un raccouci pratique. Il existe aussi -= ; *= ; /= ...
	// 	
	// 	my_prod[i] = max(NumericVector::create(20, my_prod[i]));
	// 	// juste histoire de montrer comment marche la fonction max() ici. 
	// 	// Il faut lui fournir un NumericVector donc pour ça il faut le créer avec
	// 	// NumericVector::create(x,y); (x et y étant dans scalaires)
	// 	// Ici on prendra la maximum entre 20 et la production calculée juste pour l'ex.
	// }
	// 
	// // pour la boucle while ci dessous rien de très nouveau 
	// // tu devrais tout reconnaître.
	// // la boucle ne sert à rien dans cette fonction
	// int j = 0;
	// while (j < 100) {
	// 	j++;
	// 	if(j >T1_x.size()){
	// 		break;
	// 	}
	// }
	
	ET_T2.push_back(T2_prod, "production");
	
	// Rcout << "The value of b : " << T1_inf.size() << "\n";
	// ici on ajoute le vecteur my_prod comme nouvelle (et dernière) colonne de ET_T2,
	//  grace à la fonction push_back. On peut spécifier le nom, ici "production"
	return ET_T2;
}



