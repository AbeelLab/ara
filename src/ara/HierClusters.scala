package ara

import java.io.File
import java.io.PrintWriter
import atk.compbio.tree._
import java.io.BufferedReader
import java.io.FileReader
import scala.io.Source
import scala.collection.JavaConversions._
import scala.collection.immutable.Map
import java.util.Calendar

object HierClusters {

  case class Config(
    val tree: String = null,
    val lineage: File = null,
    val out: String = null
  )

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar hier-clusters") {
      opt[String]('t', "tree") required () action { (x, c) => c.copy(tree = x) } text ("Tree file.")
      opt[File]('l', "lineage") required () action { (x, c) => c.copy(lineage = x) } text { "Macaw lineage info file." }
      opt[String]('o', "output") required () action { (x, c) => c.copy(out = x) } text ("Output prefix.")
    }

    parser.parse(args, Config()) map { config =>

      /** Configure */
      val t = new Tree(config.tree)
      val hierClusters = t.getLeaves(t.root).map(_.getName).toList
      val lin = Source.fromFile(config.lineage).getLines.filterNot { _.startsWith("#") }.map { line =>
        val arr = line.split("\t")
        (arr(1) -> arr(2))
      }.toMap.filter {
        _ match {
          case (s, lin) => hierClusters.contains(s)
        }
      } + (("reference", "LIN-4"))
      new File(config.out).mkdir()
      
      val totalTaxa = t.getLeafCount
      val linNumbers = lin.map(_._2).groupBy(identity).mapValues(_.size)
      var hc: scala.collection.immutable.Map[String, List[String]] = hierClusters.map(s => (s -> List.empty[String])).toMap

      /** Function to update hierarchical clustering map hc. */
      def updateHC[A, B](map: Map[A, List[B]], key: A, value: B) = map + ((key, map.getOrElse(key, List()) :+ value))

      /**  Function to read hierarchical clusters from a phylogenetic tree. */
      def readRecursive(node: TreeNode, lvl: Int): Unit = {
        val children = node.numberLeaves
        val child1 = node.getChild(0)
        val child2 = node.getChild(1)
        println("Parent node at level " + lvl + ": " + children + " samples")

        def getCluster(samples: List[String]): String = {
          samples.map{s => 
            val c = hc(s)
            val cName = if (c.last.contains("-")) c.last else c.filterNot(_.contains("-")).mkString(".")
            (s, cName)
          }.groupBy(_._2).keys.mkString
        }
        
        def printCluster(cluster: String, samples: List[String]): Unit = {
          val cpw = new PrintWriter(new File(config.out + "/" + cluster + "_cluster"))
          cpw.println("# Hierarchical clusters ")
          cpw.println("# " + cluster + "(N=" + samples.size + ")")
          cpw.println("# Command: java -jar ara.jar hier-clusters " + args.mkString(" "))
          cpw.println("# Compiled: " + Calendar.getInstance.getTime)
          samples.foreach(s => cpw.println(s + "\t" + cluster))
          hierClusters.diff(samples).foreach(s => cpw.println(s + "\tnot" + cluster))          
          cpw.close
        }
        
        /** Child 1 */
        if (child1.numberLeaves >= 10) {
          val c1 = t.getLeaves(child1).map(_.getName).toList
          val lineage = c1.map(lin(_)).distinct.toList
          if (child2.numberLeaves >= 10) {
            println("\tc1: " + lineage.mkString(", ") + ",  at level " + (lvl + 1) + ": " + child1.numberLeaves + " samples")
            
            val cName = if (lineage.size > 1) {
              val l = lineage.map(l => "L" + l.last).sorted.mkString("-")
              l
            } else if (linNumbers(lineage.mkString) == c1.size) {
              "L" + lineage.mkString.last
            } else {
              "1"
            }            
            c1.foreach(s => hc = updateHC(hc, s, cName))            
            
            val cluster = getCluster(c1)            
            println("\t" + cluster)
            printCluster(cluster, c1)
            
          }
          else {
            println("\tc: " + lineage.mkString(", ") + ",  at level " + (lvl + 1) + ": " + child1.numberLeaves + " samples")
          }
        } else {
          val c0 = if (child1.isLeaf()) { List(child1.getName) } else { (t.getLeaves(child1).map(_.getName)).toList }
          println("\texcluded: " + c0.mkString(", "))
        }

        /** Child 2 */
        if (child2.numberLeaves >= 10) {
          val c2 = t.getLeaves(child2).map(_.getName).toList
          val lineage = c2.map(lin(_)).distinct.toList
          if (child1.numberLeaves >= 10) {
            println("\tc2: " + lineage.mkString(", ") + ",  at level " + (lvl + 1) + ": " + child2.numberLeaves + " samples")
            
            val cName = if (lineage.size > 1) {
              val l = lineage.map(l => "L" + l.last).sorted.mkString("-")
              l
            } else if (linNumbers(lineage.mkString) == c2.size) {
              "L" + lineage.mkString.last
            } else {
              "2"
            }            
            c2.foreach(s => hc = updateHC(hc, s, cName))
            
            val cluster = getCluster(c2)            
            println("\t" + cluster)
            printCluster(cluster, c2)
            
          } else {
            println("\tc: " + lineage.mkString(", ") + ",  at level " + (lvl + 1) + ": " + child2.numberLeaves + " samples")
          }
        } else {
          val c0 = if (child2.isLeaf()) { List(child2.getName) } else { (t.getLeaves(child2).map(_.getName)).toList }
          println("\texcluded: " + c0.mkString(", "))
        }
        
        /** Print hierarchical clusters to file if both groups are defined as hierarchical clusters (size > 10).*/        
        /**if (child1.numberLeaves >= 10 && child2.numberLeaves >= 10) {
          val samples = t.getLeaves(child1).map(_.getName) ++ t.getLeaves(child2).map(_.getName).toList
          val clusters = samples.map{s => 
            val c = hc(s)
            val cName = if (c.last.contains("-")) c.last else c.filterNot(_.contains("-")).mkString(".")
            (s, cName)
          }.groupBy(_._2).mapValues(s => s.map(_._1))          
          val cpw = new PrintWriter(new File(config.out + "/" + clusters.keys.mkString("_") + "_clusters"))
          cpw.println("# Hierarchical clusters ")
          cpw.println("# " + clusters.map(_ match { case (c, ls) => (c + "(N=" + ls.size + ")")}).mkString("\n# "))
          cpw.println("# Command: java -jar ara.jar hier-clusters " + args.mkString(" "))
          cpw.println("# Compiled: " + Calendar.getInstance.getTime)
          clusters.foreach{ _ match 
            {
              case (c, ls) => ls.foreach(s => cpw.println(s + "\t" + c))
            }            
          }
          cpw.close
        }*/

        /** Recursive call */
        if (child1.numberLeaves > 19) { readRecursive(child1, lvl + 1) }
        if (child2.numberLeaves > 19) { readRecursive(child2, lvl + 1) }
      }

      /** Call function to create hierarchical clusters. */
      readRecursive(t.root, 0)

      /** Print all samples and all clusters to a file with user-defined output file name. */
      val pw = new PrintWriter(new File(config.out + "_clusters"))
      pw.println("# Hierarchical clusters")
      pw.println("# Total samples: " + totalTaxa)
      pw.println("# Command: java -jar ara.jar hier-clusters " + args.mkString(" "))
      pw.println("# Compiled: " + Calendar.getInstance.getTime)
      hc.foreach {
        _ match { case (s, ls) => pw.println(s + "\t" + ls.mkString("\t")) }
      }
      pw.close
      
      /** Print numbers of samples in each lineage. */
      println
      println("Number of samples per lineage: ")
      linNumbers.foreach(_ match {case (l, n) => println("\t" + l + " (" + n + ")")})
      
      
    }

  }

}