package ara

import java.io.File
import java.io.PrintWriter
import atk.compbio.tree._
import java.io.BufferedReader
import java.io.FileReader
import scala.io.Source
import scala.collection.JavaConversions._
import scala.collection.immutable.Map

object HierClusters {

  case class Config(
    val tree: String = null,
    val lineage: File = null,
    val out: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar hier-clusters") {
      opt[String]('t', "tree") required () action { (x, c) => c.copy(tree = x) } text ("Tree file.")
      opt[File]('l', "lineage") required () action { (x, c) => c.copy(lineage = x) } text { "Macaw lineage info file." }
      opt[File]('o', "output") required () action { (x, c) => c.copy(out = x) } text ("Output file.")
    }

    parser.parse(args, Config()) map { config =>

      /** Configure */
      val t = new Tree(config.tree)
      val hierClusters = t.getLeaves(t.root).map(_.getName).toList
      val lin = Source.fromFile(config.lineage).getLines.filterNot { _.startsWith("#") }.map { line =>
        val arr = line.split("\t")
        (arr(1) -> arr(2))
      }.toMap.filter{ _ match {
        case (s, lin) => hierClusters.contains(s)
      }} + (("reference", "LIN-4")) 
      
      val pw = new PrintWriter(config.out)

      val totalTaxa = t.getLeafCount

      
      val linNumbers = lin.map(_._2).groupBy(identity).mapValues(_.size)
      
      var hc: scala.collection.immutable.Map[String, List[String]] = hierClusters.map(s => (s -> List.empty[String])).toMap
      
      def updateHC[A, B](map: Map[A, List[B]], key: A, value: B) = map + ((key, map.getOrElse(key, List()) :+ value))
      
      def readRecursive(node: TreeNode, lvl: Int): Unit = {
        val children = node.numberLeaves
        val child1 = node.getChild(0)
        val child2 = node.getChild(1)
        println("Parent node at level " + lvl + ": " + children + " samples")
        
        /** Child 1*/
        if (child1.numberLeaves >= 10) {
          val c1 = t.getLeaves(child1).map(_.getName)
          val lineage = c1.map(lin(_)).distinct.toList
          println("\tc1: " + lineage.mkString(", ") + ",  at level " + (lvl + 1) + ": " + child1.numberLeaves + " samples")
          if (lineage.size > 1) {
            val l = lineage.sorted.mkString("/")
            c1.foreach(s => hc = updateHC(hc, s, l)) 
          }
          else if (linNumbers(lineage.mkString) == c1.size) {
            c1.foreach(s => hc = updateHC(hc, s, lineage.mkString)) 
          }
          else {
            c1.foreach(s => hc = updateHC(hc, s, "1"))            
          }          
        }
        else {
          val c0 = if (child1.isLeaf()) {List(child1.getName)} else {(t.getLeaves(child1).map(_.getName)).toList}
          println("\texcluded: " + c0.mkString(", "))
        }     
        
        /** Child 2*/
        if (child2.numberLeaves >= 10) {
          val c2 = t.getLeaves(child2).map(_.getName)
          val lineage = c2.map(lin(_)).distinct.toList
          if (lineage.size > 1) {
            val l = lineage.sorted.mkString("/")
            c2.foreach(s => hc = updateHC(hc, s, l)) 
          }
          else if (linNumbers(lineage.mkString) == c2.size) {
            c2.foreach(s => hc = updateHC(hc, s, lineage.mkString)) 
          }
          else if (child1.numberLeaves >= 10) {
            c2.foreach(s => hc = updateHC(hc, s, "2"))
            println("\tc2: " + lineage.mkString(", ") + ",  at level " + (lvl + 1) + ": " + child2.numberLeaves + " samples")
            } 
          else {
            c2.foreach(s => hc = updateHC(hc, s, "1"))
            println("\tc1: " + lineage.mkString(", ") + ",  at level " + (lvl + 1) + ": " + child2.numberLeaves + " samples")
          }
        }
        else {
          val c0 = if (child2.isLeaf()) {List(child2.getName)} else {(t.getLeaves(child2).map(_.getName)).toList}
          println("\texcluded: " + c0.mkString(", "))
        }
        
        if (child1.numberLeaves > 19) { readRecursive(child1, lvl + 1) }        
        if (child2.numberLeaves > 19) { readRecursive(child2, lvl + 1) }
        
      }

      readRecursive(t.root, 0)
      
      hc.foreach{_ match 
        { case (s, ls) => pw.println(s + "\t" + ls.mkString("\t")) }        
      }
      
      //pw.println("\n" + totalTaxa)
      pw.close
linNumbers.foreach(println)
      
    }

  }

}