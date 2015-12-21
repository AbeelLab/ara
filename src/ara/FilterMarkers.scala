package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source
import atk.util.Tool
import ara.Cluster._
import java.util.Calendar

object FilterMarkers extends Tool {
  case class Config(
    val matrix: File = null, val markerFile: File = null, output: String = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar filter-markers") {
      opt[File]("matrix") required () action { (x, c) => c.copy(matrix = x) } text ("Matrix with marker counts.")
      opt[File]("markers") required () action { (x, c) => c.copy(markerFile = x) } text ("Marker fasta file")
      opt[String]('o', "output") required () action { (x, c) => c.copy(output = x + ".txt") } text ("Output name.")
    }

    parser.parse(args, Config()) map { config =>

      /* Load matrix */
      val matrix = tLines(config.matrix)
      val matrixValues = matrix.drop(1).map(row => row.split("\t").drop(1).map(_.toInt).toList)
      println("Matrix: " + matrixValues.size + "x" + matrixValues.head.size)

      /* Load marker sequences */
      val markerSeq = tLines(config.markerFile).grouped(2).map(f => (f(0).substring(1)-> f(1))).toMap
      
      /* Read markers from first column */
      val markers = matrix.drop(1).map(_.split("\t")(0))
      println("Markers: " + markers.size)

      /* Read clusters from matrix header */
      val clusterNames = matrix.head.split("\t").drop(1).toList
      println("Clusters: " + clusterNames.size)
      println

      def rowIndex(marker: String): Int = {
        markers.indexOf(marker)
      }

      def columnIndex(cluster: String): Int = {
        clusterNames.indexOf(cluster)
      }

      def presentInCluster(marker: String): Boolean = {
        val rn = rowIndex(marker)
        val cluster = marker.split("_")(2)
        val cn = columnIndex(cluster)
        matrixValues(rn)(cn) > 0
      }

      def ancestors(cluster: String): List[String] = {
        def getParents(c: String): List[String] = {
          val ancestor = c.getAncestor
          ancestor match {
            case "MTBC" => Nil
            case x if (x.isCluster) => x :: getParents(x)
            case _ => Nil
          }          
        }
        if (cluster.isCluster && cluster != "MTBC") getParents(cluster)
        else List.empty[String]
      }

      def descendants(cluster: String): List[String] = {
        var ls = List.empty[String]
        def getChildren(c: String): Unit = {
          val children = c.children
          val c1 = children(0)
          val c2 = children(1)
          if (c1.isCluster) { ls = c1 :: ls; getChildren(c1) }
          if (c2.isCluster) { ls = c2 :: ls; getChildren(c2) }
        }
        getChildren(cluster)
        ls.sorted
      }

      def presentInUnrelatedCluster(marker: String): Boolean = {
        val rn = rowIndex(marker)
        val cluster = marker.split("_")(2)
        val related = ancestors(cluster) ++ descendants(cluster)
        val unrelated = clusterNames.filterNot { c => related.contains(c) }
        val count = unrelated.map{c =>
          val cn = columnIndex(c)
          matrixValues(rn)(cn)
        }.foldLeft(0)(_ + _)
        count > 0        
      }

      var countA = 0
      var countB = 0
      
      val filteredMarkers = markers.filter{m => 
        val rn = rowIndex(m)
        if (rn % 1000 == 0 && rn != 0) println(rn + " markers checked..")
        val a = presentInCluster(m) 
        val b = !presentInUnrelatedCluster(m)     
        if (a) countA = countA + 1
        if (b) countB = countB + 1
        a && b
      }
      
      println("Markers detected in specific cluster: " + countA)
      println("Markers not detected in unrelated clusters: " + countB)
      println("Filtered markers: " + filteredMarkers.size)
      
      val pw = new PrintWriter(config.output)
      pw.println("## Filtered markers after validation")
      pw.println("## Command: " + args.mkString(" "))
      pw.println("## Date: " + Calendar.getInstance.getTime)
      pw.println("#")
      filteredMarkers.sortBy{m => m.split("_")(0).drop(1).dropRight(1).toInt}.foreach{m => pw.println(">" + m + "/n" + markerSeq(m))}
      pw.close

    }

  }
}