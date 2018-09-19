/***************************************************************************/
/* gauss.hpp                                                               */
/*                                                                         */
/*   Diese Datei definiert Konstanten fuer die GaussQuadratur              */
/*                                                                         */
/*   Gaus Auswertungspunkte auf dem Intervall [0 1]                        */
/*                                                                         */
/*   2^$NUM entspricht der Anzahl der Auswertungspunkte                    */
/*                                                                         */
/*   Anzahl der Punkte:    GAUSS_SIZE[$NUM]                                */
/*   Auswertungstellen:    GAUSS_NODES[$NUM][$I].n  ($I aktuelle Position) */
/*   Gewichte:             GAUSS_NODES[$NUM][$I].w  ($I aktuelle Position) */
/*                                                                         */
/*                                                                         */
/*   Beispiel:                                                             */
/*                                                                         */
/*   > double sol;                                                         */
/*   > for(int i=0;i<GAUSS_SIZE[3];++i)                                    */
/*   >    sol += sin(GAUSS_NODES[3][i].n) * GAUSS_NODES[3][i].w;           */
/*                                                                         */
/***************************************************************************/
/* Author: Peter Schaefer                             schaeferpm@gmail.com */
/* Version: 1.0  (2012)                                                    */
/***************************************************************************/

#ifndef GAUSS_NODES_
#define GAUSS_NODES_

typedef struct _gauss {
	double n;
	double w;
} gauss;

double GAUSS_SIZE[] = { 1, 2, 4, 8, 16, 32 };

gauss GAUSS_NODES[][32] = {
{
	{ 0.5, 1 }},
{
	{ 0.211324865405187078959, 0.499999999999999888978 },
	{ 0.788675134594812865529, 0.499999999999999888978 }},
{
	{ 0.0694318442029737137311, 0.173927422568727063634 },
	{ 0.330009478207571815833, 0.32607257743127315841 },
	{ 0.669990521792428128656, 0.326072577431272991877 },
	{ 0.93056815579702634178, 0.173927422568727008123 }},
{
	{ 0.0198550717512319119251, 0.0506142681451882195387 },
	{ 0.101666761293186636017, 0.111190517226687282659 },
	{ 0.237233795041835615613, 0.15685332293894366229 },
	{ 0.408282678752175054449, 0.181341891689181078373 },
	{ 0.591717321247825056574, 0.18134189168918116164 },
	{ 0.762766204958164384387, 0.156853322938943606779 },
	{ 0.898333238706813475005, 0.111190517226687282659 },
	{ 0.980144928248767977053, 0.0506142681451882542332 }},
{
	{ 0.00529953250417525278948, 0.0135762297058772823943 },
	{ 0.0277124884633836443548, 0.0311267619693235096656 },
	{ 0.0671843988060840668908, 0.0475792558412463095774 },
	{ 0.12229779582249833414, 0.0623144856277671119194 },
	{ 0.191061877798677837159, 0.0747979944082881875733 },
	{ 0.270991611171386259649, 0.0845782596975013928331 },
	{ 0.359198224610370597798, 0.0913017075224623053664 },
	{ 0.452493745081181397705, 0.0947253052275341400623 },
	{ 0.547506254918818657806, 0.0947253052275344453736 },
	{ 0.640801775389629457713, 0.0913017075224619306661 },
	{ 0.729008388828613518307, 0.0845782596975014205887 },
	{ 0.808938122201321996307, 0.0747979944082887010515 },
	{ 0.877702204177501554838, 0.0623144856277674172307 },
	{ 0.932815601193915933109, 0.0475792558412462263107 },
	{ 0.972287511536616300134, 0.0311267619693240300827 },
	{ 0.994700467495824858233, 0.0135762297058771956582 }},
{
	{ 0.00136806907525910403933, 0.00350930500473500784839 },
	{ 0.00719424422736580915227, 0.00813719736545263742922 },
	{ 0.0176188722062468050567, 0.0126960326546312531754 },
	{ 0.0325469620311301111037, 0.0171369314565104867432 },
	{ 0.0518394221169737878796, 0.0214179490111136711095 },
	{ 0.0753161931337150702959, 0.0254990296311881116387 },
	{ 0.102758102016028862735, 0.0293420467392677027096 },
	{ 0.133908940629855199855, 0.0329111113881807790249 },
	{ 0.168477866534892439798, 0.0361728970544242522944 },
	{ 0.206142121379618625809, 0.0390969478935349543103 },
	{ 0.246550045533885375804, 0.041655962113473471442 },
	{ 0.289324361934682361408, 0.0438260465022020373471 },
	{ 0.334065698858936110938, 0.0455869393478821535726 },
	{ 0.380356318873931620317, 0.0469221995404019084908 },
	{ 0.427764019208601742328, 0.0478193600396371459871 },
	{ 0.475846167156130650522, 0.0482700442573640656208 },
	{ 0.524153832843869071922, 0.0482700442573634549981 },
	{ 0.572235980791398368694, 0.0478193600396370696592 },
	{ 0.619643681126068490705, 0.046922199540402220741 },
	{ 0.66593430114106377804, 0.0455869393478817580556 },
	{ 0.710675638065317860637, 0.0438260465022019540804 },
	{ 0.753449954466114735219, 0.041655962113473277153 },
	{ 0.793857878620381374191, 0.0390969478935357661609 },
	{ 0.831522133465107504691, 0.0361728970544242176 },
	{ 0.866091059370144966678, 0.032911111388180640247 },
	{ 0.897241897983971248287, 0.029342046739267681893 },
	{ 0.924683806866285151749, 0.0254990296311881359248 },
	{ 0.94816057788302621212, 0.021417949011113403962 },
	{ 0.967453037968869833385, 0.0171369314565111112436 },
	{ 0.982381127793753194943, 0.0126960326546313780754 },
	{ 0.992805755772633968803, 0.00813719736545309192677 },
	{ 0.998631930924740673916, 0.00350930500473496925079 }}};

#endif