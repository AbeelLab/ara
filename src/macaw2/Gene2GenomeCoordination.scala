package macaw2

import scala.io.Source
import java.io.File

/**
 * @author Arlin
 */
object Gene2GenomeCoordination {

  val usage = "scala Gene2GenomeCoordination.scala [list] [gff-file] [fasta-file]"

  def main(args: Array[String]) {

    if (args.size != 3) {
      println(usage)
    } else {

            
      val ls = Source.fromFile(new File(args(0))).getLines.filterNot(_.startsWith("#")).map { line =>
        val arr = line.split("\t")
        val nChange = arr(4).split("/")
        val oldRefBases = nChange.take(nChange.size/2).map(_.charAt(0)).mkString("/")
        (arr(8), arr(3), oldRefBases)
      }.toList

      class Gene(val tag: String, val start: Int, val end: Int, val dir: String) {
        override def toString(): String = "Gene(" + tag + ", " + start + ", " + end + ", " + dir + ")"
      }

      object Gene {
        def unapply(s: String): Option[(String, Int, Int, String)] = {
          val arr = s.split("\t")
          val source = arr(2)
          val start = arr(3).toInt
          val end = arr(4).toInt
          val direction = arr(6)
          val locus_tag = arr(8).split(";").filter(_.startsWith("locus_tag")).mkString.split("=")(1)
          Some(locus_tag.substring(1, locus_tag.size - 1), start, end, direction)
        }
      }

      val gff = Source.fromFile(new File(args(1))).getLines.filter(line => line.split("\t")(2) != "CDS" && line.split("\t")(2) != "source").map(line => line match {
        case Gene(tag, start, end, dir) => new Gene(tag, start, end, dir)
      }).toList
      
      val genes = gff.filter(g => ls.map(_._1).contains(g.tag)).map(g => (g.tag -> g)).toMap
      
      val ref = Source.fromFile(new File(args(2))).getLines.filterNot(_.startsWith(">")).mkString 
      
      val newLs = ls.map{_ match {
        case (locus_tag, geneCoor, oldRefBases) => {
          val gene = genes(locus_tag)
          geneCoor match {
            case s if (s != "-") => {
              val arr = s.split("/")
              val gC = arr.map(_.toInt)
              val chrCoor = gC.map{_ match {
                case c if c < 0 => {
                  if (gene.dir == "+" ) 
                    if (locus_tag == "RVBD_6018") gene.start + c - 10 
                    else if (locus_tag == "RVBD_6019") gene.start + c - 2
                    else (gene.start + c)
                  else (gene.end - c)  
                }
                case c if c > 0 => {
                  if (gene.dir == "+" ) 
                    if (locus_tag == "RVBD_6018") gene.start + c - 11
                    else if (locus_tag == "RVBD_6019") gene.start + c - 3
                    else (gene.start + c - 1)
                  else (gene.end - c + 1)  
                }                
              }}
              val chrBases = chrCoor.map(c => ref.charAt(c - 1)).mkString("/")
              (locus_tag, geneCoor, chrCoor.mkString("/"), chrBases, chrBases == oldRefBases)              
            }
            case _ => (locus_tag, geneCoor, "-", "-", "-")            
          }          
        }
      }}
      
      newLs.foreach(m => println(m._3 + "\t" + m._2 + "\t" + m._4 + "\t" + m._1 + "\t" + m._5))
      
    }
  }
}