package macaw2

import java.io.File
import scala.io.Source
import java.io.PrintWriter


object DrugResistances {
  
  case class Config(
    val snpTyperOutput: String = null, 
    val result: String = null
  )
  
  class Marker(val mutInfo: String, val count: Int, val isPresent: Boolean){
    override def toString(): String = mutInfo + "\t" + count + "\t" + (if (isPresent) "present" else "absent")
    val arr = mutInfo.split("_")
    val coordinate = arr(0).drop(1).dropRight(1).toInt
    val markerType = arr(1)
    val drugs = arr(2)
  }

  object Marker {
    def unapply(s: String): Option[(String, Int, Boolean)] = {
      val line = s.split("\t")
      val mutInfo = line(0)
      val count = line(1).toInt
      val presence = if (line(3) == "present") true else false
      Some(mutInfo, count, presence)
    }
  }
  
  def main(args: Array[String]){
    
    val parser = new scopt.OptionParser[Config]("java -jar Ara.jar interpret-DR"){
      opt[String]('m', "markers") required() action { (x, c) => c.copy(snpTyperOutput = x) } text ("Output file of MacawSNPTyper.")
      opt[String]('o', "output") required() action { (x, c) => c.copy(result = x) } text ("Output name for the file with results.")
    }
    
    class findInList(val ls: List[Marker]){
      val presentMarkers = ls.filter(_.isPresent)
      def hasResistance = {
        if (presentMarkers.exists(_.markerType == "resistance")) true
        else false
      }
      def hasSusceptibility = {
        if (presentMarkers.exists(_.markerType == "susceptibility")) true
        else false
      }
    }
    
    implicit def checkMarkers(ls: List[Marker]) = new findInList(ls)
    
    parser.parse(args, Config()) map { config =>
      val snpTypes = new File(config.snpTyperOutput)
      val outFile = new File(config.result)
      
      val pw = new PrintWriter(outFile)
      pw.println("# Results drug resistance")
      
      val markers = Source.fromFile(snpTypes).getLines.filterNot(_.startsWith("#")).toList.dropRight(2).map {_ match {
          case Marker(mutInfo, count, isPresent) => new Marker(mutInfo, count, isPresent)
      }}.groupBy(_.coordinate)//.map(_._2).toList.sortBy(_(1).coordinate)
      
      //markers.foreach(pw.println)
      pw.println(markers.size + " genome positions to detect that are associated with drug susceptibility or resistance.\n")
      
            
      val resistances = markers.filterNot(_._2.hasSusceptibility).filter(_._2.hasResistance).map(_._2.filter(_.isPresent))
      pw.println(resistances.size + " genome position(s) detected marking only drug resistance:\n")
      resistances.foreach(l => pw.println("\t" + l.mkString("\n\t") + "\n"))
      
      val susceptibilities = markers.filterNot(_._2.hasResistance).filter(_._2.hasSusceptibility).map(_._2.filter(_.isPresent))
      //susceptibilities.foreach(pw.println)
      pw.println(susceptibilities.size + " genome position(s) detected marking only drug susceptibility.\n")
      
      val nonDetected = markers.filterNot(_._2.hasResistance).filterNot(_._2.hasSusceptibility)
      pw.println(nonDetected.size + " non-detected SNP positions with possible resistance:\n")
      
      def ranges(ls: List[Int], range: Int): List[List[Int]] = ls match {
        case head :: tail => {
          val (b, a) = ls.span(p => p - head < 11)
          b:: ranges(a, 10)
        }
        case Nil => Nil
      }
      
      val within10bp = ranges(nonDetected.map(_._1).toList.sorted, 10).map(_.flatMap(nonDetected(_)))
      within10bp.foreach(l => pw.println("\t" + l.mkString("\n\t") + "\n"))
      
      
      
      within10bp.map(_.filter(_.markerType == "susceptibility")).foreach{ pList =>
        val drugs = pList.flatMap(_.drugs.split("_").toList).distinct
        val coordinates = pList.map(_.coordinate)
        if (coordinates.size > 1) pw.println("\tPossible resistance to " + drugs.mkString(", ") + " caused by one or more SNPs at genome coordinations " + coordinates.mkString(", ") + ".")
        else pw.println("\tPossible resistance to " + drugs.mkString(", ") + " caused by a SNP at genome coordination " + coordinates.mkString + " that is not included in the drug resistance list.")
      }
      pw.println
      
      val bothDetected = markers.filter(_._2.hasResistance).filter(_._2.hasSusceptibility)
      pw.println(bothDetected.size + " genome positions with both resistance and susceptibility detected:\n")
      bothDetected.foreach(l => pw.println("\t" + l._2.mkString("\n\t") + "\n"))
      
      pw.close
      
      
    }
    
  }
}