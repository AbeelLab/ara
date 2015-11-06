package ara


object Cluster {
  class Cluster(s: String){
    
    val refCluster = "L4.2.2.1.1.2.2.1.1.1.1.1"
    
    def hasReference: Boolean = {        
       (s.contains("L4") || s.equals("L4") || s.equals(refCluster.dropRight(10)) || s.equals(refCluster.dropRight(8)) || s.equals(refCluster.dropRight(6)) || 
              s.equals(refCluster.dropRight(4)) || s.equals(refCluster.dropRight(2)) || s.equals(refCluster))
      }
  }
  implicit def seqtoBool(s: String) = new Cluster(s)
}