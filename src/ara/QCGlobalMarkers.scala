package ara

import java.io.File
import java.io.PrintWriter
import scala.io.Source
import ara.Cluster._
import atk.util.Tool
import Array._

object QCGlobalMarkers extends MTBCclusters with Tool {

  case class Config(
    val markerFile: File = null,
    val output: String = null,
    val directory: File = null,
    val clusterDir: File = null,
    val trnSet: File = null)

  def main(args: Array[String]) {

    val parser = new scopt.OptionParser[Config]("java -jar ara.jar markers-qc") {
      opt[File]("markers") action { (x, c) => c.copy(markerFile = x) } text ("File containing marker sequences. This file has to be a multi-fasta file with the headers indicating the name of the markers.") //, { v: String => config.spacerFile = v })
      opt[File]('o', "output") required () action { (x, c) => c.copy(output = x + "_qc.txt") } text ("Output name.")
      opt[File]('d', "data-dir") required () action { (x, c) => c.copy(directory = x) } text ("Data directory.")
      opt[File]('c', "cluster-dir") required () action { (x, c) => c.copy(clusterDir = x) } text ("Directory with cluster files.")
      opt[File]('t', "trn-set") required () action { (x, c) => c.copy(trnSet = x) } text ("VCF path file")
    }

    parser.parse(args, Config()) map { config =>

      def listAraFiles(f: Any): List[File] = f match {
        case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listAraFiles(_))
        case f: File if (f.isFile() && f.getName.endsWith(".ara") && !f.getName.endsWith("interpret.ara")) => List(f)
        case _ => Nil
      }

      def listClusterFiles(f: Any): List[File] = f match {
        case f: File if (f.isDirectory()) => f.listFiles().toList.flatMap(listClusterFiles(_))
        case f: File if (f.isFile() && f.getName.endsWith("_cluster")) => List(f)
        case _ => Nil
      }

      /** Elapsed time function */
      def time[R](block: => R): R = {
        val t0 = System.currentTimeMillis()
        val result = block // call-by-name
        val t1 = System.currentTimeMillis()
        println("Elapsed time: " + (t1 - t0) + "ms")
        result
      }

      time {
        /* Load markers */
        val lines = if (config.markerFile != null) tLines(config.markerFile).toList else scala.io.Source.fromInputStream(MacawSNPtyper.getClass().getResourceAsStream("/hierclusters_global_bluejaymarkers")).getLines().filterNot(f => f.startsWith("#") || f.trim.size == 0).toList
        val markers = lines.grouped(2).map(f => (f(0).substring(1))).toList.sortBy(_.split("_")(2))
        println("Markers: " + markers.size)
        /* Load cluster files */
        val clusterNames = markers.map(_.split("_")(2)).distinct.sorted
        val clusterMap = listClusterFiles(config.clusterDir).filter(f => clusterNames.contains(f.getName.dropRight(8))).map { c =>
          val strains = tLines(c).filterNot(_.contains("not")).map(_.split("\t")(0))
          strains.map(s => (s, c.getName.split("_")(0)))
        }.flatten.groupBy(_._1).mapValues(_.map(_._2))
        println("Clusters with  markersets: " + clusterNames.size)
        /* Load SNP-typer files */
        val araFiles = listAraFiles(config.directory).map{f => 
          val study = f.getParentFile.getParentFile.getName
          val sampleID = f.getParentFile.getName 
          (study, sampleID) -> f
        }.toMap//.filter(f => trnset.contains(f.getParentFile.getName))
        val trnset = tLines(config.trnSet).map{path => 
          val arr = path.split("/")
          val study = arr(7)
          val sampleID = arr(8)
          (study, sampleID)
        }.filter(araFiles.keysIterator.toList.contains(_)).map(araFiles(_))
        println("Trainingset: " + araFiles.size)
        println

        /* Prepare matrix */
        val rowhead = "$$" :: clusterNames
        var matrix = ofDim[Int](markers.size, clusterNames.size)

        def rowIndex(marker: String): Int = {
          markers.indexOf(marker)
        }

        def columnIndex(cluster: String): Int = {
          clusterNames.indexOf(cluster)
        }

        def readAraFile(file: File): Unit = {
          val sampleID = file.getParentFile.getName
          print(sampleID + ", ")
          val clusters = clusterMap(sampleID)
          val presentMarkers = tLines(file).dropRight(2).map { line =>
            line match {
              case Marker(m, c, p) => new ClusterMarker(m, c, p)
            }
          }.filter(_.isPresent)
          println(presentMarkers.size + " present markers")
          presentMarkers.foreach { marker =>
            for (cn <- clusters) {
              val prevCount = matrix(rowIndex(marker.mutInfo))(columnIndex(cn))
              matrix(rowIndex(marker.mutInfo))(columnIndex(cn)) = prevCount + 1
            }
          }
        }

        trnset.foreach { readAraFile(_) }

        def printMatrix(m: Array[Array[Int]]): Unit = {
          val pw = new PrintWriter(config.output)
          pw.println(rowhead.mkString("\t"))
          for (rn <- 0 until markers.size) {
            pw.print(markers(rn))
            for (cn <- 0 until clusterNames.size) {
              pw.print("\t" + m(rn)(cn))
            }
            pw.println
          }
          pw.close
        }
        printMatrix(matrix)
      }

    }
  }
}